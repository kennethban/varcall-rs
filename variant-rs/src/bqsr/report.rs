//! The recalibration report: finalizing the count tables (`BaseRecalibrationEngine.finalizeData`,
//! `QuantizationInfo`) and GATK's `GATKReport` v1.1 text format, byte for byte.

use super::covariates::{context_from_key, cycle_from_key, ReadGroups, CONTEXT_DIM, CYCLE_DIM, QUAL_DIM};
use super::datum::{quantize, Datum};
use super::recal::Counts;
use hc_genotype::javamath::format_fixed;
use rayon::prelude::*;
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub enum Val {
    Str(String),
    Long(i64),
    Double(f64),
}

impl Val {
    fn cmp_java(&self, o: &Val) -> Ordering {
        match (self, o) {
            (Val::Str(a), Val::Str(b)) => a.cmp(b),
            (Val::Long(a), Val::Long(b)) => a.cmp(b),
            (Val::Double(a), Val::Double(b)) => crate::jmath::java_cmp(*a, *b),
            _ => self.to_plain().cmp(&o.to_plain()),
        }
    }
    fn to_plain(&self) -> String {
        match self {
            Val::Str(s) => s.clone(),
            Val::Long(v) => v.to_string(),
            Val::Double(v) => format!("{v}"),
        }
    }
}

pub struct Table {
    pub name: String,
    pub desc: String,
    /// (column name, format: "%s", "%d", "%.4f", "%.2f", or "" for an untyped column)
    pub cols: Vec<(String, String)>,
    pub rows: Vec<Vec<Val>>,
}

/// `NumberUtils.isCreatable(s.trim())` for the values that can appear in a recalibration report.
fn is_creatable(s: &str) -> bool {
    let s = s.trim();
    let b = s.as_bytes();
    if b.is_empty() {
        return false;
    }
    let mut i = 0;
    if b[0] == b'-' || b[0] == b'+' {
        i = 1;
    }
    if i >= b.len() {
        return false;
    }
    let mut digits = 0;
    let mut dot = false;
    let mut exp = false;
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_digit() {
            digits += 1;
        } else if c == b'.' && !dot && !exp {
            dot = true;
        } else if (c == b'e' || c == b'E') && digits > 0 && !exp {
            exp = true;
            if i + 1 < b.len() && (b[i + 1] == b'-' || b[i + 1] == b'+') {
                i += 1;
            }
            if i + 1 >= b.len() {
                return false;
            }
        } else if i == b.len() - 1 && digits > 0 && matches!(c, b'd' | b'D' | b'f' | b'F') {
            return true;
        } else if i == b.len() - 1 && digits > 0 && !dot && !exp && matches!(c, b'l' | b'L') {
            return true;
        } else {
            return false;
        }
        i += 1;
    }
    digits > 0
}

fn is_right_align(v: &str) -> bool {
    matches!(v, "null" | "NA" | "Infinity" | "-Infinity" | "NaN") || is_creatable(v)
}

fn format_val(v: &Val, fmt: &str) -> String {
    match v {
        Val::Str(s) => s.clone(),
        Val::Long(x) => x.to_string(),
        Val::Double(x) => {
            if fmt.is_empty() {
                format_fixed(*x, 8)
            } else if !x.is_finite() {
                if x.is_nan() { "NaN".into() } else if *x > 0.0 { "Infinity".into() } else { "-Infinity".into() }
            } else {
                let prec: usize = fmt.trim_start_matches("%.").trim_end_matches('f').parse().expect("format");
                format_fixed(*x, prec)
            }
        }
    }
}

impl Table {
    /// `GATKReportTable.write` with SORT_BY_COLUMN.
    pub fn write(&mut self, out: &mut String) {
        let rows = &mut self.rows;
        rows.sort_by(|a, b| {
            for (x, y) in a.iter().zip(b.iter()) {
                let c = x.cmp_java(y);
                if c != Ordering::Equal {
                    return c;
                }
            }
            Ordering::Equal
        });
        let nc = self.cols.len();
        let cells: Vec<Vec<String>> = rows.iter().map(|r| r.iter().enumerate().map(|(i, v)| format_val(v, &self.cols[i].1)).collect()).collect();
        let mut width: Vec<usize> = self.cols.iter().map(|c| c.0.len()).collect();
        let mut right = vec![true; nc];
        for r in &cells {
            for (i, s) in r.iter().enumerate() {
                if !s.is_empty() {
                    width[i] = width[i].max(s.len());
                    if right[i] && !is_right_align(s) {
                        right[i] = false;
                    }
                }
            }
        }
        out.push_str(&format!("#:GATKTable:{}:{}", nc, rows.len()));
        for c in &self.cols {
            out.push(':');
            out.push_str(if c.1.is_empty() { "%s" } else { &c.1 });
        }
        out.push_str(":;\n");
        out.push_str(&format!("#:GATKTable:{}:{}\n", self.name, self.desc));
        for (i, c) in self.cols.iter().enumerate() {
            if i > 0 {
                out.push_str("  ");
            }
            out.push_str(&format!("{:<w$}", c.0, w = width[i]));
        }
        out.push('\n');
        for r in &cells {
            for (i, s) in r.iter().enumerate() {
                if i > 0 {
                    out.push_str("  ");
                }
                if right[i] {
                    out.push_str(&format!("{:>w$}", s, w = width[i]));
                } else {
                    out.push_str(&format!("{:<w$}", s, w = width[i]));
                }
            }
            out.push('\n');
        }
        out.push('\n');
    }
}

/// Finalized tables (after collapse and rounding), keyed like GATK's NestedIntegerArrays.
pub struct Final {
    pub rg: Vec<Option<Datum>>,                       // [rg]
    pub qual: Vec<Option<Datum>>,                     // [rg * QUAL_DIM + q]
    pub context: Vec<(usize, usize, usize, Datum)>,   // (rg, q, key, datum) in key order
    pub cycle: Vec<(usize, usize, usize, Datum)>,
    pub quantized: Vec<u8>,
    pub qual_hist: Vec<i64>,
}

/// `finalizeData()` + `QuantizationInfo(tables, quantizingLevels)`.
pub fn finalize(c: &Counts, quantizing_levels: usize) -> Final {
    let n_rg = c.n_rg;
    let mut qual: Vec<Option<Datum>> = c
        .qual
        .iter()
        .enumerate()
        .map(|(i, x)| if x[0] > 0 { Some(Datum::from_counts(x[0], x[1], (i % QUAL_DIM) as u8)) } else { None })
        .collect();
    // collapseQualityScoreTableToReadGroupTable: leaves of the qual table in key order
    let mut rg: Vec<Option<Datum>> = vec![None; n_rg];
    for (i, d) in qual.iter().enumerate() {
        if let Some(d) = d {
            let r = i / QUAL_DIM;
            match &mut rg[r] {
                None => rg[r] = Some(*d),
                Some(x) => x.combine(d),
            }
        }
    }
    let collect = |v: &Vec<[u64; 2]>, dim: usize| -> Vec<(usize, usize, usize, Datum)> {
        v.par_iter()
            .enumerate()
            .filter(|(_, x)| x[0] > 0)
            .map(|(i, x)| {
                let key = i % dim;
                let rq = i / dim;
                let q = rq % QUAL_DIM;
                (rq / QUAL_DIM, q, key, Datum::from_counts(x[0], x[1], q as u8))
            })
            .collect()
    };
    let mut context = collect(&c.context, CONTEXT_DIM);
    let mut cycle = collect(&c.cycle, CYCLE_DIM);
    // roundTableValues
    for d in rg.iter_mut().chain(qual.iter_mut()).flatten() {
        d.round();
    }
    for t in [&mut context, &mut cycle] {
        for x in t.iter_mut() {
            x.3.round();
        }
    }
    // QuantizationInfo: histogram of the qual table's empirical qualities (prior = reported)
    let mut hist = vec![0i64; QUAL_DIM];
    for d in qual.iter().flatten() {
        let e = d.empirical_quality(d.reported);
        hist[e as usize] += d.obs;
    }
    let quantized = quantize(&hist, quantizing_levels, 6);
    Final { rg, qual, context, cycle, quantized, qual_hist: hist }
}

fn ev() -> Val {
    Val::Str("M".into())
}

/// The whole report text (`RecalUtils.outputRecalibrationReport`).
pub fn report_text(f: &Final, rgs: &ReadGroups, quantizing_levels: usize) -> String {
    report_text_with(f, &rgs.names, arguments_table(quantizing_levels))
}

/// BaseRecalibrator's Arguments table (default arguments).
fn arguments_table(quantizing_levels: usize) -> Table {
    let args: [(&str, &str); 17] = [
        ("covariate", "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate"),
        ("no_standard_covs", "false"),
        ("run_without_dbsnp", "false"),
        ("solid_recal_mode", "SET_Q_ZERO"),
        ("solid_nocall_strategy", "THROW_EXCEPTION"),
        ("mismatches_context_size", "2"),
        ("indels_context_size", "3"),
        ("mismatches_default_quality", "-1"),
        ("deletions_default_quality", "45"),
        ("insertions_default_quality", "45"),
        ("maximum_cycle_value", "500"),
        ("low_quality_tail", "2"),
        ("default_platform", "null"),
        ("force_platform", "null"),
        ("quantizing_levels", ""),
        ("recalibration_report", "null"),
        ("binary_tag_name", "null"),
    ];
    let ql = quantizing_levels.to_string();
    Table {
        name: "Arguments".into(),
        desc: "Recalibration argument collection values used in this run".into(),
        cols: vec![("Argument".into(), "%s".into()), ("Value".into(), "".into())],
        rows: args.iter().map(|(a, v)| vec![Val::Str(a.to_string()), Val::Str(if *a == "quantizing_levels" { ql.clone() } else { v.to_string() })]).collect(),
    }
}

/// The report text with the given Arguments table and read-group names (`RecalUtils.createRecalibrationGATKReport`).
pub fn report_text_with(f: &Final, names: &[String], mut args: Table) -> String {
    let mut out = String::from("#:GATKReport.v1.1:5\n");
    args.write(&mut out);
    // Quantized
    let mut t = Table {
        name: "Quantized".into(),
        desc: "Quality quantization map".into(),
        cols: vec![("QualityScore".into(), "%d".into()), ("Count".into(), "%d".into()), ("QuantizedScore".into(), "%d".into())],
        rows: (0..QUAL_DIM).map(|q| vec![Val::Long(q as i64), Val::Long(f.qual_hist[q]), Val::Long(f.quantized[q] as i64)]).collect(),
    };
    t.write(&mut out);
    // RecalTable0: read groups
    let mut t = Table {
        name: "RecalTable0".into(),
        desc: "".into(),
        cols: vec![
            ("ReadGroup".into(), "%s".into()),
            ("EventType".into(), "%s".into()),
            ("EmpiricalQuality".into(), "%.4f".into()),
            ("EstimatedQReported".into(), "%.4f".into()),
            ("Observations".into(), "%d".into()),
            ("Errors".into(), "%.2f".into()),
        ],
        rows: vec![],
    };
    for (r, d) in f.rg.iter().enumerate() {
        if let Some(d) = d {
            t.rows.push(vec![
                Val::Str(names[r].clone()),
                ev(),
                Val::Double(d.empirical_quality(d.reported) as f64),
                Val::Double(d.reported),
                Val::Long(d.obs),
                Val::Double(d.num_mismatches()),
            ]);
        }
    }
    t.write(&mut out);
    // RecalTable1: quality scores
    let mut t = Table {
        name: "RecalTable1".into(),
        desc: "".into(),
        cols: vec![
            ("ReadGroup".into(), "%s".into()),
            ("QualityScore".into(), "%d".into()),
            ("EventType".into(), "%s".into()),
            ("EmpiricalQuality".into(), "%.4f".into()),
            ("Observations".into(), "%d".into()),
            ("Errors".into(), "%.2f".into()),
        ],
        rows: vec![],
    };
    for (i, d) in f.qual.iter().enumerate() {
        if let Some(d) = d {
            t.rows.push(vec![
                Val::Str(names[i / QUAL_DIM].clone()),
                Val::Long((i % QUAL_DIM) as i64),
                ev(),
                Val::Double(d.empirical_quality(d.reported) as f64),
                Val::Long(d.obs),
                Val::Double(d.num_mismatches()),
            ]);
        }
    }
    t.write(&mut out);
    // RecalTable2: context and cycle
    let mut t = Table {
        name: "RecalTable2".into(),
        desc: "".into(),
        cols: vec![
            ("ReadGroup".into(), "%s".into()),
            ("QualityScore".into(), "%d".into()),
            ("CovariateValue".into(), "%s".into()),
            ("CovariateName".into(), "%s".into()),
            ("EventType".into(), "%s".into()),
            ("EmpiricalQuality".into(), "%.4f".into()),
            ("Observations".into(), "%d".into()),
            ("Errors".into(), "%.2f".into()),
        ],
        rows: vec![],
    };
    let rows: Vec<Vec<Val>> = f
        .context
        .par_iter()
        .map(|x| (x, true))
        .chain(f.cycle.par_iter().map(|x| (x, false)))
        .map(|(&(r, q, key, d), is_ctx)| {
            vec![
                Val::Str(names[r].clone()),
                Val::Long(q as i64),
                Val::Str(if is_ctx { context_from_key(key) } else { cycle_from_key(key).to_string() }),
                Val::Str(if is_ctx { "Context".into() } else { "Cycle".into() }),
                ev(),
                Val::Double(d.empirical_quality(d.reported) as f64),
                Val::Long(d.obs),
                Val::Double(d.num_mismatches()),
            ]
        })
        .collect();
    t.rows = rows;
    t.write(&mut out);
    out
}
