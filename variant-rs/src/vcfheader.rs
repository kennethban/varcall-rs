//! htsjdk's repair of standard VCF header lines, which its VCF reader applies to every header it reads; tools that
//! pass an input's header lines through (MergeVcfs, ApplyVQSR) apply it too.

/// htsjdk's standard header lines (`VCFStandardHeaderLines`): (FORMAT?, ID, Number, Type, Description).
const STANDARD_LINES: &[(bool, &str, &str, &str, &str)] = &[
    (true, "GT", "1", "String", "Genotype"),
    (true, "GQ", "1", "Integer", "Genotype Quality"),
    (true, "DP", "1", "Integer", "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"),
    (true, "PL", "G", "Integer", "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"),
    (true, "AD", "R", "Integer", "Allelic depths for the ref and alt alleles in the order listed"),
    (true, "FT", ".", "String", "Genotype-level filter"),
    (true, "PS", "1", "Integer", "Phasing set (typically the position of the first variant in the set)"),
    (true, "PQ", "1", "Float", "Read-backed phasing quality"),
    (false, "END", "1", "Integer", "Stop position of the interval"),
    (false, "DB", "0", "Flag", "dbSNP Membership"),
    (false, "DP", "1", "Integer", "Approximate read depth; some reads may have been filtered"),
    (false, "SB", "1", "Float", "Strand Bias"),
    (false, "AF", "A", "Float", "Allele Frequency, for each ALT allele, in the same order as listed"),
    (false, "AC", "A", "Integer", "Allele count in genotypes, for each ALT allele, in the same order as listed"),
    (false, "AN", "1", "Integer", "Total number of alleles in called genotypes"),
    (false, "MQ0", "1", "Integer", "Total Mapping Quality Zero Reads"),
    (false, "MQ", "1", "Float", "RMS Mapping Quality"),
    (false, "SOMATIC", "0", "Flag", "Somatic event"),
];

/// The `key=value` items of a `##KEY=<...>` line's body (after `<`; a final `>` dropped), commas inside quotes kept.
pub fn compound_fields(inner: &str) -> Vec<String> {
    let inner = inner.strip_suffix('>').unwrap_or(inner);
    let (mut fields, mut cur, mut quoted, mut escaped) = (Vec::new(), String::new(), false, false);
    for c in inner.chars() {
        match c {
            _ if escaped => escaped = false,
            '\\' => escaped = true,
            '"' => quoted = !quoted,
            ',' if !quoted => {
                fields.push(std::mem::take(&mut cur));
                continue;
            }
            _ => {}
        }
        cur.push(c);
    }
    fields.push(cur);
    fields
}

/// (key, ID, Number, Type) of an `##INFO` / `##FORMAT` / `##FILTER` / `##contig` style line (`None` for other lines).
pub fn compound(line: &str) -> Option<(String, String, Option<String>, Option<String>, Vec<String>)> {
    let body = line.strip_prefix("##")?;
    let (key, rest) = body.split_once("=<")?;
    let fields = compound_fields(rest);
    let get = |k: &str| fields.iter().find_map(|f| f.strip_prefix(k).and_then(|r| r.strip_prefix('='))).map(|s| s.to_string());
    Some((key.to_string(), get("ID")?, get("Number"), get("Type"), fields))
}

/// `VCFStandardHeaderLines.repairStandardHeaderLines`, which htsjdk's VCF reader applies to every header it reads: an
/// INFO / FORMAT line for a standard field whose Number (count type, or count when fixed) or Type differs from the
/// standard is replaced by the standard line (e.g. `AD` with `Number=.` becomes `Number=R`). Other lines are kept.
pub fn repair_standard_line(line: String) -> String {
    let (format, inner) = if let Some(r) = line.strip_prefix("##FORMAT=<") {
        (true, r)
    } else if let Some(r) = line.strip_prefix("##INFO=<") {
        (false, r)
    } else {
        return line;
    };
    let fields = compound_fields(inner);
    let get = |k: &str| fields.iter().find_map(|f| f.strip_prefix(k).and_then(|r| r.strip_prefix('=')));
    let (Some(id), Some(number), Some(typ)) = (get("ID"), get("Number"), get("Type")) else { return line };
    let Some(&(_, _, std_number, std_type, desc)) = STANDARD_LINES.iter().find(|s| s.0 == format && s.1 == id) else { return line };
    // count type: A, R, G, "." (unbounded) or a fixed integer; a fixed count is compared only when both are fixed
    let bad_count = match (number.parse::<i64>(), std_number.parse::<i64>()) {
        (Ok(n), Ok(s)) => n != s,
        (Err(_), Err(_)) => number != std_number,
        _ => true,
    };
    if bad_count || typ != std_type {
        format!("##{}=<ID={id},Number={std_number},Type={std_type},Description=\"{desc}\">", if format { "FORMAT" } else { "INFO" })
    } else {
        line
    }
}

/// htsjdk's `samplesWereAlreadySorted`. When a VCF's sample names are not in sorted order, htsjdk's reader decodes
/// every genotype and its writer re-encodes them (FORMAT keys recomputed, values reformatted), which the tools here
/// do not reproduce: they write the genotype columns as read, so they refuse such inputs.
pub fn samples_sorted(samples: &[String]) -> bool {
    samples.windows(2).all(|w| w[0].encode_utf16().cmp(w[1].encode_utf16()) != std::cmp::Ordering::Greater)
}

/// Refuses an input whose samples are not sorted (see `samples_sorted`).
pub fn require_sorted_samples(samples: &[String], path: &str) {
    if !samples_sorted(samples) {
        crate::gatk_cli::refuse(&format!(
            "{path}: the sample columns are not in sorted order; htsjdk then re-encodes every genotype, which is not supported"
        ));
    }
}

/// The output header as htsjdk's VCF writer prints it: `##fileformat=VCFv4.2`, then the input's lines and `extra`
/// (exact duplicates once) sorted as text, the contig lines as one block in their order where the first sorts, and
/// our command line where GATK's `##GATKCommandLine` would sort; then the #CHROM line with `samples`.
pub fn write_header(input_meta: &[String], extra: Vec<String>, samples: &[String]) -> String {
    let mut lines: Vec<String> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for l in input_meta.iter().cloned().chain(extra) {
        if l.starts_with("##fileformat=") {
            continue;
        }
        if seen.insert(l.clone()) {
            lines.push(l);
        }
    }
    let key = |l: &String| -> String {
        if l.starts_with("##contig=<") {
            "##contig=<".into()
        } else {
            l.strip_prefix("##variant-rsCommandLine").map(|r| format!("##GATKCommandLine{r}")).unwrap_or_else(|| l.clone())
        }
    };
    lines.sort_by_key(key);
    let mut out = String::from("##fileformat=VCFv4.2\n");
    for l in lines {
        out.push_str(&l);
        out.push('\n');
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !samples.is_empty() {
        out.push_str("\tFORMAT");
        for s in samples {
            out.push('\t');
            out.push_str(s);
        }
    }
    out.push('\n');
    out
}

/// htsjdk's `escapeQuotes` for a quoted header value: a '"' not preceded by a backslash gets one (as
/// `replaceAll("([^\\\\])\"", ...)`, so a leading quote, or the second of two, is left alone).
fn escape_quotes(v: &str) -> String {
    let c: Vec<char> = v.chars().collect();
    let mut out = String::with_capacity(v.len() + 4);
    let mut i = 0;
    while i < c.len() {
        if i + 1 < c.len() && c[i] != '\\' && c[i + 1] == '"' {
            out.push(c[i]);
            out.push_str("\\\"");
            i += 2;
        } else {
            out.push(c[i]);
            i += 1;
        }
    }
    out
}

/// A `##FILTER` line as `VCFFilterHeaderLine` writes it.
pub fn filter_line(id: &str, description: &str) -> String {
    format!("##FILTER=<ID={id},Description=\"{}\">", escape_quotes(description))
}

/// This run's command line, written where GATK writes `##GATKCommandLine` (see `write_header`).
pub fn command_line(tool: &str) -> String {
    let cmdline: Vec<String> = std::env::args().collect();
    format!(
        "##variant-rsCommandLine=<ID={tool},CommandLine=\"{}\",Version=\"{}\",Date=\"{}\">",
        cmdline.join(" ").replace('"', "\\\""),
        env!("CARGO_PKG_VERSION"),
        crate::now_string()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_lines() {
        let ad = "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths, here unbounded\">";
        assert_eq!(repair_standard_line(ad.into()), "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
        // only the description differs: kept
        let dp = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">";
        assert_eq!(repair_standard_line(dp.into()), dp);
        // the type differs
        let mq = "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"x, y\">";
        assert_eq!(repair_standard_line(mq.into()), "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">");
        // not a standard field, or the INFO/FORMAT namesake of one
        let qd = "##INFO=<ID=QD,Number=.,Type=String,Description=\"x\">";
        assert_eq!(repair_standard_line(qd.into()), qd);
        let gq = "##INFO=<ID=GQ,Number=.,Type=String,Description=\"x\">";
        assert_eq!(repair_standard_line(gq.into()), gq);
    }

    #[test]
    fn quotes() {
        assert_eq!(escape_quotes(r#"TYPE == "SNP""#), r#"TYPE == \"SNP\""#);
        assert_eq!(escape_quotes(r#""a" x"#), r#""a\" x"#);
        assert_eq!(escape_quotes(r#"a "" b"#), r#"a \"" b"#);
        assert_eq!(escape_quotes(r#"a \" b"#), r#"a \" b"#);
    }
}
