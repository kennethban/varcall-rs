//! Truth-sensitivity tranches (`TrancheManager.findTranches`, `Tranche`, `TruthSensitivityTranche`): finding them
//! on the lod-sorted data, the tranches file VariantRecalibrator writes, and reading it back for ApplyVQSR.

use super::{Datum, Mode};
use hc_genotype::javamath::format_fixed;

pub struct Tranche {
    pub target: f64,
    pub min_vqslod: f64,
    pub num_known: i64,
    pub known_titv: f64,
    pub num_novel: i64,
    pub novel_titv: f64,
    pub accessible_truth_sites: i64,
    pub calls_at_truth_sites: i64,
    pub model: String,
    pub filter_name: String,
}

impl Tranche {
    fn truth_sensitivity(&self) -> f64 {
        if self.accessible_truth_sites > 0 { self.calls_at_truth_sites as f64 / (1.0 * self.accessible_truth_sites as f64) } else { 0.0 }
    }

    /// `getTrancheString`: `%.2f,%d,%d,%.4f,%.4f,%.4f,VQSRTranche%s%.2fto%.2f,%s,%d,%d,%.4f%n`
    fn line(&self, prev: Option<&Tranche>) -> String {
        format!(
            "{},{},{},{},{},{},VQSRTranche{}{}to{},{},{},{},{}\n",
            format_fixed(self.target, 2),
            self.num_known,
            self.num_novel,
            format_fixed(self.known_titv, 4),
            format_fixed(self.novel_titv, 4),
            format_fixed(self.min_vqslod, 4),
            self.model,
            format_fixed(prev.map(|p| p.target).unwrap_or(0.0), 2),
            format_fixed(self.target, 2),
            self.model,
            self.accessible_truth_sites,
            self.calls_at_truth_sites,
            format_fixed(self.truth_sensitivity(), 4)
        )
    }
}

fn count_calls_at_truth(data: &[Datum], min_lod: f64) -> i64 {
    data.iter().filter(|d| d.at_truth && d.lod >= min_lod).count() as i64
}

/// `Tranche.trancheOfVariants` + `TruthSensitivityTranche`
fn tranche_of_variants(data: &[Datum], min_i: usize, ts: f64, mode: Mode) -> Tranche {
    let (mut num_known, mut num_novel) = (0i64, 0i64);
    let (mut known_ti, mut known_tv, mut novel_ti, mut novel_tv) = (0i64, 0i64, 0i64, 0i64);
    let min_lod = data[min_i].lod;
    for d in data {
        if d.lod >= min_lod {
            if d.is_known {
                num_known += 1;
                if d.is_snp {
                    if d.is_transition { known_ti += 1 } else { known_tv += 1 }
                }
            } else {
                num_novel += 1;
                if d.is_snp {
                    if d.is_transition { novel_ti += 1 } else { novel_tv += 1 }
                }
            }
        }
    }
    Tranche {
        target: ts,
        min_vqslod: min_lod,
        num_known,
        known_titv: known_ti as f64 / (1.0 * known_tv as f64).max(1.0),
        num_novel,
        novel_titv: novel_ti as f64 / (1.0 * novel_tv as f64).max(1.0),
        accessible_truth_sites: count_calls_at_truth(data, f64::NEG_INFINITY),
        calls_at_truth_sites: count_calls_at_truth(data, min_lod),
        model: mode.name().to_string(),
        filter_name: "anonymous".to_string(),
    }
}

/// Sorts `data` by lod (as `findTranches` does, stably), finds the tranches and returns the tranches file text.
pub fn tranches_text(data: &mut [Datum], thresholds: &[f64], mode: Mode) -> Result<String, String> {
    let n_true_sites = count_calls_at_truth(data, f64::NEG_INFINITY);
    data.sort_by(|a, b| hc_assembly::jmath::java_cmp(a.lod, b.lod));
    // running truth sensitivity from the top
    let mut running = vec![0.0; data.len()];
    let mut n_called = 0i64;
    for i in (0..data.len()).rev() {
        n_called += data[i].at_truth as i64;
        running[i] = 1.0 - n_called as f64 / (1.0 * n_true_sites as f64);
    }
    let mut tranches: Vec<Tranche> = Vec::new();
    for &t in thresholds {
        let metric_threshold = 1.0 - t / 100.0;
        match (0..data.len()).find(|&i| running[i] >= metric_threshold) {
            Some(i) => tranches.push(tranche_of_variants(data, i, t, mode)),
            None => {
                if tranches.is_empty() {
                    return Err(format!(
                        "Couldn't find any tranche containing variants with a TruthSensitivity > {}. Are you sure the truth files contain unfiltered variants which overlap the input data?",
                        format_fixed(metric_threshold, 2)
                    ));
                }
                break;
            }
        }
    }
    // Tranche.tranchesString: sorted by calls at truth sites (stably), each line knowing the previous one
    tranches.sort_by(|a, b| a.calls_at_truth_sites.cmp(&b.calls_at_truth_sites));
    let mut s = String::from("# Variant quality score tranches file\n# Version number 5\n");
    s.push_str("targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n");
    let mut prev: Option<&Tranche> = None;
    for t in &tranches {
        s.push_str(&t.line(prev));
        prev = Some(t);
    }
    Ok(s)
}

/// `TruthSensitivityTranche.readTranches`, sorted by target truth sensitivity.
pub fn read_tranches(path: &str) -> Result<Vec<Tranche>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut header: Option<Vec<String>> = None;
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let vals: Vec<&str> = line.split(',').collect();
        let Some(h) = &header else {
            if vals.len() != 11 {
                return Err(format!("{path}: Expected 11 elements in header line {line}"));
            }
            header = Some(vals.iter().map(|s| s.to_string()).collect());
            continue;
        };
        if vals.len() != h.len() {
            return Err(format!("{path}: Line had too few/many fields.  Header = {} vals {}. The line was: {line}", h.len(), vals.len()));
        }
        let get = |k: &str| h.iter().position(|x| x == k).map(|i| vals[i]);
        let dbl = |k: &str, def: Option<f64>| -> Result<f64, String> {
            match get(k) {
                Some(v) => super::vcfio::parse_vcf_double(v).ok_or(format!("Malformed tranches file. Invalid value for key {k}")),
                None => def.ok_or(format!("Malformed tranches file.  Missing required key {k}")),
            }
        };
        let int = |k: &str, def: Option<i64>| -> Result<i64, String> {
            match get(k) {
                Some(v) => v.parse::<i32>().map(|x| x as i64).map_err(|_| format!("Malformed tranches file. Invalid value for key {k}")),
                None => def.ok_or(format!("Malformed tranches file.  Missing required key {k}")),
            }
        };
        out.push(Tranche {
            target: dbl("targetTruthSensitivity", None)?,
            min_vqslod: dbl("minVQSLod", None)?,
            num_known: int("numKnown", Some(-1))?,
            known_titv: dbl("knownTiTv", Some(-1.0))?,
            num_novel: int("numNovel", None)?,
            novel_titv: dbl("novelTiTv", None)?,
            accessible_truth_sites: int("accessibleTruthSites", Some(-1))?,
            calls_at_truth_sites: int("callsAtTruthSites", Some(-1))?,
            model: get("model").unwrap_or("").to_string(),
            filter_name: get("filterName").unwrap_or("").to_string(),
        });
    }
    out.sort_by(|a, b| hc_assembly::jmath::java_cmp(a.target, b.target));
    Ok(out)
}
