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
    let inner = inner.strip_suffix('>').unwrap_or(inner);
    // key=value pairs, commas inside quotes kept
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
}
