//! The two gather steps of a scattered Mutect2 run (nf-core/sarek runs both once per tumor):
//! - `MergeMutectStats`: the `callable` statistic of the scatters' `.stats` tables, summed;
//! - `GatherPileupSummaries`: the scatters' GetPileupSummaries tables, the empty ones dropped, the rest ordered by
//!   their first site (contig index in the dictionary, then position; a stable sort) and concatenated.

use super::contam::{read_pileups, write_pileups};
use super::fmc_math::java_double_to_string;

fn fail(msg: String) -> ! {
    eprintln!("A USER ERROR has occurred: {msg}");
    std::process::exit(2);
}

/// `DoubleStream.sum()` of JDK 17 (GATK's runtime): Kahan summation plus a simple sum, `computeFinalSum` adding the
/// compensation term (JDK 18 subtracts it).
fn java17_double_stream_sum(v: &[f64]) -> f64 {
    let (mut sum, mut comp, mut simple) = (0.0f64, 0.0f64, 0.0f64);
    for &d in v {
        let tmp = d - comp;
        let velvel = sum + tmp;
        comp = (velvel - sum) - tmp;
        sum = velvel;
        simple += d;
    }
    let tmp = sum + comp;
    if tmp.is_nan() && simple.is_infinite() {
        simple
    } else {
        tmp
    }
}

const MMS_USAGE: &str = "usage: variant-rs mergemutectstats --stats a.stats [--stats b.stats ...] -O merged.stats";

/// `MergeMutectStats`: the input set keeps the first of repeated paths (`LinkedHashSet<File>`); every statistic must
/// be `callable` (the only aggregation GATK defines).
pub fn run_merge_stats(argv: &[String]) {
    let (mut stats, mut out) = (Vec::<String>::new(), None::<String>);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{MMS_USAGE}"));
        match argv[i].as_str() {
            "--stats" => {
                let p = v();
                if !stats.contains(&p) {
                    stats.push(p);
                }
            }
            "-O" | "--output" => out = Some(v()),
            a => panic!("unknown argument {a}\n{MMS_USAGE}"),
        }
        i += 2;
    }
    let out = out.unwrap_or_else(|| panic!("{MMS_USAGE}"));
    let mut callable: Vec<f64> = Vec::new();
    for p in &stats {
        let text = std::fs::read_to_string(p).unwrap_or_else(|e| fail(format!("Encountered an IO exception while reading from {p}: {e}")));
        let mut lines = text.lines().filter(|l| !l.starts_with('#'));
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let (ks, vs) = match (header.iter().position(|&h| h == "statistic"), header.iter().position(|&h| h == "value")) {
            (Some(k), Some(v)) => (k, v),
            _ => fail(format!("{p}: not a Mutect2 stats table (columns statistic, value)")),
        };
        for l in lines.filter(|l| !l.is_empty()) {
            let f: Vec<&str> = l.split('\t').collect();
            if f[ks] != "callable" {
                fail(format!("aggregations list missing key {}", f[ks]));
            }
            callable.push(f[vs].parse().unwrap_or_else(|_| fail(format!("{p}: bad value {}", f[vs]))));
        }
    }
    let mut text = String::from("statistic\tvalue\n");
    if !callable.is_empty() {
        text.push_str(&format!("callable\t{}\n", java_double_to_string(java17_double_stream_sum(&callable))));
    }
    std::fs::write(&out, text).unwrap_or_else(|e| fail(format!("Encountered an IO exception while writing to {out}: {e}")));
}

const GPS_GATHER_USAGE: &str = "usage: variant-rs gatherpileupsummaries -I a.pileups.table [-I b.pileups.table ...] -O out.table --sequence-dictionary ref.dict";

/// `GatherPileupSummaries`.
pub fn run_gather_pileups(argv: &[String]) {
    let (mut inputs, mut out, mut dict) = (Vec::<String>::new(), None::<String>, None::<String>);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{GPS_GATHER_USAGE}"));
        match argv[i].as_str() {
            "-I" | "--I" => inputs.push(v()),
            "-O" | "--O" => out = Some(v()),
            "--sequence-dictionary" => dict = Some(v()),
            a => panic!("unknown argument {a}\n{GPS_GATHER_USAGE}"),
        }
        i += 2;
    }
    let (out, dict) = (out.unwrap_or_else(|| panic!("{GPS_GATHER_USAGE}")), dict.unwrap_or_else(|| panic!("{GPS_GATHER_USAGE}")));
    // ReferenceUtils.loadFastaDictionary: the @SQ names in order
    let dict_text = std::fs::read_to_string(&dict).unwrap_or_else(|e| fail(format!("{dict}: {e}")));
    let contigs: Vec<&str> = dict_text.lines().filter(|l| l.starts_with("@SQ")).filter_map(|l| l.split('\t').find_map(|f| f.strip_prefix("SN:"))).collect();
    // PileupSummaryComparator: List.indexOf (-1 for a contig not in the dictionary), then the start
    let index_of = |c: &str| contigs.iter().position(|&x| x == c).map_or(-1, |k| k as i64);
    let mut tables: Vec<(String, Vec<super::contam::Site>)> = inputs.iter().map(|p| read_pileups(p)).filter(|(_, s)| !s.is_empty()).collect();
    tables.sort_by_key(|(_, s)| (index_of(&s[0].contig), s[0].pos)); // stable, as Collections.sort
    if tables.is_empty() {
        std::fs::write(&out, "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n").unwrap_or_else(|e| fail(format!("{out}: {e}")));
        return;
    }
    let sample = tables[0].0.clone();
    if let Some((s, _)) = tables.iter().find(|(s, _)| *s != sample) {
        fail(format!("Combining PileupSummaryTables from different samples is not supported. Got samples {sample} and {s}"));
    }
    let sites: Vec<super::contam::Site> = tables.into_iter().flat_map(|(_, s)| s).collect();
    write_pileups(&out, &sample, &sites);
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn jdk17_sum() {
        assert_eq!(java17_double_stream_sum(&[1.0, 2.0, 3.0]), 6.0);
        assert_eq!(java17_double_stream_sum(&[]), 0.0);
        // values from JDK 17's DoubleStream.of(...).sum() (jshell)
        assert_eq!(java17_double_stream_sum(&[0.1, 0.2, 0.3]), 0.6);
        assert_eq!(java17_double_stream_sum(&[1e16, 1.0, -1e16]), 0.0);
        assert_eq!(java17_double_stream_sum(&[0.1; 10]), 1.0);
    }
}
