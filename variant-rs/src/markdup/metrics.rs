//! Picard `DuplicationMetrics` and htsjdk `MetricsFile` output (FormatUtil number formatting).

use std::collections::BTreeMap;

#[derive(Default, Clone, Debug)]
pub struct LibMetrics {
    pub unpaired_reads_examined: i64,
    pub read_pairs_examined: i64, // counted per read until finalized (then / 2)
    pub secondary_or_supplementary_rds: i64,
    pub unmapped_reads: i64,
    pub unpaired_read_duplicates: i64,
    pub read_pair_duplicates: i64, // per read until finalized
    pub read_pair_optical_duplicates: i64,
}

impl LibMetrics {
    pub fn add(&mut self, o: &LibMetrics) {
        self.unpaired_reads_examined += o.unpaired_reads_examined;
        self.read_pairs_examined += o.read_pairs_examined;
        self.secondary_or_supplementary_rds += o.secondary_or_supplementary_rds;
        self.unmapped_reads += o.unmapped_reads;
        self.unpaired_read_duplicates += o.unpaired_read_duplicates;
        self.read_pair_duplicates += o.read_pair_duplicates;
        self.read_pair_optical_duplicates += o.read_pair_optical_duplicates;
    }

    /// `DuplicationMetrics.addReadToLibraryMetrics` (+ `addDuplicateReadToMetrics` when `dup`).
    pub fn add_read(&mut self, flag: u16, dup: bool) {
        let unmapped = flag & 0x4 != 0;
        let sec_supp = flag & 0x900 != 0;
        let unpaired = flag & 0x1 == 0 || flag & 0x8 != 0;
        if unmapped {
            self.unmapped_reads += 1;
        } else if sec_supp {
            self.secondary_or_supplementary_rds += 1;
        } else if unpaired {
            self.unpaired_reads_examined += 1;
        } else {
            self.read_pairs_examined += 1;
        }
        if dup && !sec_supp && !unmapped {
            if unpaired {
                self.unpaired_read_duplicates += 1;
            } else {
                self.read_pair_duplicates += 1;
            }
        }
    }
}

fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + crate::jmath::exp(-n / x)
}

/// `DuplicationMetrics.estimateLibrarySize`.
pub fn estimate_library_size(read_pairs: i64, unique_read_pairs: i64) -> Option<i64> {
    let dups = read_pairs - unique_read_pairs;
    if read_pairs > 0 && dups > 0 {
        let (c, n) = (unique_read_pairs as f64, read_pairs as f64);
        let mut m = 1.0f64;
        let mut mm = 100.0f64;
        assert!(!(unique_read_pairs >= read_pairs || f(m * c, c, n) < 0.0), "Invalid values for pairs and unique pairs: {read_pairs}, {unique_read_pairs}");
        while f(mm * c, c, n) > 0.0 {
            mm *= 10.0;
        }
        for _ in 0..40 {
            let r = (m + mm) / 2.0;
            let u = f(r * c, c, n);
            if u == 0.0 {
                break;
            } else if u > 0.0 {
                m = r;
            } else if u < 0.0 {
                mm = r;
            }
        }
        Some((c * (m + mm) / 2.0) as i64)
    } else {
        None
    }
}

/// htsjdk `FormatUtil.format(double)`: `NumberFormat` without grouping, at most 6 fraction digits,
/// HALF_DOWN on the exact binary value, trailing zeros dropped; NaN/Infinity print as "?".
pub fn format_double(x: f64) -> String {
    if !x.is_finite() {
        return if x.is_nan() { "?".into() } else if x > 0.0 { "?".into() } else { "-?".into() };
    }
    let neg = x < 0.0;
    let exact = format!("{:.80}", x.abs()); // exact decimal expansion (doubles have <= 1074 digits; 80 suffices here)
    let (ip, fp) = exact.split_once('.').unwrap();
    let mut digits: Vec<u8> = ip.bytes().chain(fp.bytes().take(6)).collect();
    let rest: &str = if fp.len() > 6 { &fp[6..] } else { "" };
    // HALF_DOWN: round up only if the remainder is strictly greater than half
    let first = rest.bytes().next().unwrap_or(b'0');
    let round_up = first > b'5' || (first == b'5' && rest.bytes().skip(1).any(|c| c != b'0'));
    if round_up {
        let mut i = digits.len();
        loop {
            if i == 0 {
                digits.insert(0, b'1');
                break;
            }
            i -= 1;
            if digits[i] == b'9' {
                digits[i] = b'0';
            } else {
                digits[i] += 1;
                break;
            }
        }
    }
    let n_int = digits.len() - 6;
    let int_part: String = String::from_utf8(digits[..n_int].to_vec()).unwrap().trim_start_matches('0').to_string();
    let int_part = if int_part.is_empty() { "0".to_string() } else { int_part };
    let frac: String = String::from_utf8(digits[n_int..].to_vec()).unwrap().trim_end_matches('0').to_string();
    let mut s = String::new();
    if neg && (int_part != "0" || !frac.is_empty()) {
        s.push('-');
    }
    s.push_str(&int_part);
    if !frac.is_empty() {
        s.push('.');
        s.push_str(&frac);
    }
    s
}

/// Java `Double.toString` for the integer-valued histogram bins (< 1e7).
fn java_double_key(k: i64) -> String {
    format!("{k}.0")
}

/// The metrics file: one row per library (sorted by name), ROI histogram when there is one library,
/// and the duplicate-set-size histograms.
pub fn write(
    path: &str,
    cmdline: &str,
    libs: &BTreeMap<String, LibMetrics>,
    all_sets: &BTreeMap<i64, i64>,
    optical_sets: &BTreeMap<i64, i64>,
    non_optical_sets: &BTreeMap<i64, i64>,
) {
    let mut out = String::new();
    out.push_str("## htsjdk.samtools.metrics.StringHeader\n");
    out.push_str(&format!("# {cmdline}\n"));
    out.push_str("## htsjdk.samtools.metrics.StringHeader\n");
    out.push_str(&format!("# Started on: {}\n\n", std::process::Command::new("date").output().map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string()).unwrap_or_default()));
    out.push_str("## METRICS CLASS\tpicard.sam.DuplicationMetrics\n");
    out.push_str("LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n");
    let mut finals = Vec::new();
    for (name, m) in libs {
        let mut m = m.clone();
        m.read_pairs_examined /= 2;
        m.read_pair_duplicates /= 2;
        let els = estimate_library_size(m.read_pairs_examined - m.read_pair_optical_duplicates, m.read_pairs_examined - m.read_pair_duplicates);
        let pct = if m.unpaired_reads_examined + m.read_pairs_examined != 0 {
            (m.unpaired_read_duplicates + m.read_pair_duplicates * 2) as f64 / (m.unpaired_reads_examined + m.read_pairs_examined * 2) as f64
        } else {
            0.0
        };
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            name,
            m.unpaired_reads_examined,
            m.read_pairs_examined,
            m.secondary_or_supplementary_rds,
            m.unmapped_reads,
            m.unpaired_read_duplicates,
            m.read_pair_duplicates,
            m.read_pair_optical_duplicates,
            format_double(pct),
            els.map(|v| v.to_string()).unwrap_or_default()
        ));
        finals.push((m, els));
    }
    out.push('\n');
    // histograms: ROI (single library only), then all_sets, optical_sets, non_optical_sets; empty ones skipped
    let mut cols: Vec<(&str, BTreeMap<i64, f64>)> = Vec::new();
    if finals.len() == 1 {
        if let (m, Some(els)) = &finals[0] {
            let unique = m.read_pairs_examined - m.read_pair_duplicates;
            let mut h = BTreeMap::new();
            for x in 1..=100 {
                let v = *els as f64 * (1.0 - crate::jmath::exp(-(x as f64 * m.read_pairs_examined as f64) / *els as f64)) / unique as f64;
                h.insert(x, v);
            }
            cols.push(("CoverageMult", h));
        }
    }
    for (label, h) in [("all_sets", all_sets), ("optical_sets", optical_sets), ("non_optical_sets", non_optical_sets)] {
        if !h.is_empty() {
            cols.push((label, h.iter().map(|(&k, &v)| (k, v as f64)).collect()));
        }
    }
    if !cols.is_empty() {
        // htsjdk writes the first histogram's bin label: the ROI histogram's "BIN", else the duplicate-set
        // histograms' "set_size" (LibraryIdGenerator: a library without duplicate pairs has no ROI histogram)
        let bin_label = if cols[0].0 == "CoverageMult" { "BIN" } else { "set_size" };
        out.push_str("## HISTOGRAM\tjava.lang.Double\n");
        out.push_str(bin_label);
        for (l, _) in &cols {
            out.push('\t');
            out.push_str(l);
        }
        out.push('\n');
        let keys: std::collections::BTreeSet<i64> = cols.iter().flat_map(|(_, h)| h.keys().copied()).collect();
        for k in keys {
            out.push_str(&java_double_key(k));
            for (_, h) in &cols {
                out.push('\t');
                out.push_str(&format_double(*h.get(&k).unwrap_or(&0.0)));
            }
            out.push('\n');
        }
    }
    // MetricsFile.write: a newline after the histogram section, also when there is no histogram
    out.push('\n');
    std::fs::write(path, out).unwrap_or_else(|e| panic!("{path}: {e}"));
}
