//! GatherBQSRReports (`variant-rs gatherbqsrreports`): `RecalibrationReport.gatherReports`. Every input report is
//! parsed back into `RecalDatum`s (observations, the printed 2-decimal errors, and the reported quality: the
//! `EstimatedQReported` column in the read-group table, the quality score elsewhere), keyed by the sorted union of
//! all inputs' read groups; tables are combined left to right (`RecalUtils.combineTables`: a key missing so far takes
//! the other datum, else `RecalDatum.combine`), the quantization is recomputed from the combined quality-score table,
//! and the report is printed with the first input's Arguments table. No rounding after combining.

use super::covariates::{key_from_context, key_from_cycle, QUAL_DIM};
use super::datum::{quantize, Datum, MULTIPLIER};
use super::report::{report_text_with, Final, Table, Val};
use std::collections::{BTreeSet, HashMap};

const USAGE: &str = "usage: variant-rs gatherbqsrreports -I recal1.table [-I recal2.table ...] -O gathered.table";

/// One parsed `GATKReportTable`: column (name, format) and rows as text.
struct RawTable {
    cols: Vec<(String, String)>,
    desc: String,
    rows: Vec<Vec<String>>,
}

/// `GATKReport(File)` for the tables of a recalibration report.
fn parse_report(path: &str) -> HashMap<String, RawTable> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut lines = text.lines();
    let first = lines.next().unwrap_or("");
    assert!(first.starts_with("#:GATKReport.v1.1:"), "{path}: not a GATK report (v1.1)");
    let mut tables = HashMap::new();
    while let Some(l) = lines.next() {
        let Some(spec) = l.strip_prefix("#:GATKTable:") else { continue };
        // #:GATKTable:ncols:nrows:fmt...:;  then  #:GATKTable:name:description
        let f: Vec<&str> = spec.trim_end_matches(":;").split(':').collect();
        let (nc, nr): (usize, usize) = (f[0].parse().unwrap(), f[1].parse().unwrap());
        let fmts: Vec<String> = f[2..2 + nc].iter().map(|s| s.to_string()).collect();
        let title = lines.next().and_then(|l| l.strip_prefix("#:GATKTable:")).unwrap_or_else(|| panic!("{path}: table without a name"));
        let (name, desc) = title.split_once(':').unwrap_or((title, ""));
        let names: Vec<String> = lines.next().unwrap().split_whitespace().map(|s| s.to_string()).collect();
        assert_eq!(names.len(), nc, "{path}: table {name}: column names");
        let rows: Vec<Vec<String>> = (0..nr)
            .map(|_| {
                let r: Vec<String> = lines.next().unwrap().split_whitespace().map(|s| s.to_string()).collect();
                assert_eq!(r.len(), nc, "{path}: table {name}: row");
                r
            })
            .collect();
        tables.insert(name.to_string(), RawTable { cols: names.into_iter().zip(fmts).collect(), desc: desc.to_string(), rows });
    }
    tables
}

impl RawTable {
    fn col(&self, name: &str) -> usize {
        self.cols.iter().position(|c| c.0 == name).unwrap_or_else(|| panic!("column {name} missing"))
    }
}

/// Table keys: (table, read group, quality, covariate value, covariate name).
type Key = (u8, String, i32, String, String);

/// `getRecalDatum`: `new RecalDatum(obs, errors, 1)` then `setReportedQuality`.
fn datum(t: &RawTable, r: &[String], reported: f64) -> Datum {
    let obs: f64 = r[t.col("Observations")].parse().unwrap();
    let err: f64 = r[t.col("Errors")].parse().unwrap();
    Datum { obs: obs as i64, mism: err * MULTIPLIER, reported }
}

fn report_datums(path: &str, t: &HashMap<String, RawTable>) -> Vec<(Key, Datum)> {
    let mut out = Vec::new();
    let get = |n: &str| t.get(n).unwrap_or_else(|| panic!("{path}: table {n} missing"));
    let event_ok = |tb: &RawTable, r: &[String]| {
        let e = &r[tb.col("EventType")];
        assert!(e == "M", "{path}: event type {e}: variant-rs gathers base substitution tables only (BaseRecalibrator's default)");
    };
    let t0 = get("RecalTable0");
    for r in &t0.rows {
        event_ok(t0, r);
        let q: f64 = r[t0.col("EstimatedQReported")].parse().unwrap();
        out.push(((0, r[t0.col("ReadGroup")].clone(), 0, String::new(), String::new()), datum(t0, r, q)));
    }
    let t1 = get("RecalTable1");
    for r in &t1.rows {
        event_ok(t1, r);
        let q: i32 = r[t1.col("QualityScore")].parse().unwrap();
        out.push(((1, r[t1.col("ReadGroup")].clone(), q, String::new(), String::new()), datum(t1, r, q as f64)));
    }
    let t2 = get("RecalTable2");
    for r in &t2.rows {
        event_ok(t2, r);
        let q: i32 = r[t2.col("QualityScore")].parse().unwrap();
        let key = (2, r[t2.col("ReadGroup")].clone(), q, r[t2.col("CovariateValue")].clone(), r[t2.col("CovariateName")].clone());
        out.push((key, datum(t2, r, q as f64)));
    }
    out
}

pub fn run(argv: &[String]) {
    let (mut inputs, mut output) = (Vec::<String>::new(), None::<String>);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{USAGE}"));
        match argv[i].as_str() {
            "-I" | "--input" => inputs.push(v()),
            "-O" | "--output" => output = Some(v()),
            a => panic!("unknown argument {a}\n{USAGE}"),
        }
        i += 2;
    }
    let output = output.unwrap_or_else(|| panic!("{USAGE}"));
    assert!(!inputs.is_empty(), "Cannot gather an empty list of inputs");
    let reports: Vec<HashMap<String, RawTable>> = inputs.iter().map(|p| parse_report(p)).collect();
    // all read groups (TreeSet): the read-group covariate's keys
    let names: Vec<String> = reports
        .iter()
        .flat_map(|t| t["RecalTable0"].rows.iter().map(|r| r[t["RecalTable0"].col("ReadGroup")].clone()).collect::<Vec<_>>())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let rg_index: HashMap<&str, usize> = names.iter().enumerate().map(|(k, n)| (n.as_str(), k)).collect();
    // combine, left to right
    let mut combined: HashMap<Key, Datum> = HashMap::new();
    for (p, t) in inputs.iter().zip(&reports) {
        for (k, d) in report_datums(p, t) {
            match combined.get_mut(&k) {
                Some(x) => x.combine(&d),
                None => {
                    combined.insert(k, d);
                }
            }
        }
    }
    let first = &reports[0]["Arguments"];
    let levels: usize = first
        .rows
        .iter()
        .find(|r| r[0] == "quantizing_levels")
        .map(|r| r[1].parse().expect("quantizing_levels"))
        .unwrap_or(16);
    let n = names.len();
    let mut f = Final { rg: vec![None; n], qual: vec![None; n * QUAL_DIM], context: Vec::new(), cycle: Vec::new(), quantized: Vec::new(), qual_hist: Vec::new() };
    for ((tb, rg, q, val, cov), d) in combined {
        let r = rg_index[rg.as_str()];
        let q = q as usize;
        match tb {
            0 => f.rg[r] = Some(d),
            1 => f.qual[r * QUAL_DIM + q] = Some(d),
            _ => match cov.as_str() {
                "Context" => f.context.push((r, q, key_from_context(val.as_bytes()) as usize, d)),
                "Cycle" => f.cycle.push((r, q, key_from_cycle(val.parse().unwrap()), d)),
                c => panic!("covariate {c}: not a standard covariate"),
            },
        }
    }
    // QuantizationInfo(recalibrationTables, quantizingLevels)
    let mut hist = vec![0i64; QUAL_DIM];
    for d in f.qual.iter().flatten() {
        hist[d.empirical_quality(d.reported) as usize] += d.obs;
    }
    f.quantized = quantize(&hist, levels, 6);
    f.qual_hist = hist;
    let args = Table {
        name: "Arguments".into(),
        desc: first.desc.clone(),
        cols: first.cols.clone(),
        rows: first.rows.iter().map(|r| r.iter().map(|v| Val::Str(v.clone())).collect()).collect(),
    };
    let text = report_text_with(&f, &names, args);
    std::fs::write(&output, text).unwrap_or_else(|e| panic!("{output}: {e}"));
    eprintln!("gatherbqsrreports: {} reports, {} read groups -> {output}", inputs.len(), n);
}
