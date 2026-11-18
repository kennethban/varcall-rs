//! VCF data lines exactly as htsjdk's `VCFEncoder` writes them for HaplotypeCaller calls.

use crate::javamath::format_qual;

/// One output record. `info` values are already formatted; keys are sorted (htsjdk writes the
/// INFO map through a `TreeMap`).
#[derive(Clone, Debug)]
pub struct VcfRecord {
    pub contig: String,
    pub pos: i64,
    /// ID column (`--dbsnp` rsIDs); None = "."
    pub id: Option<String>,
    /// reference allele first
    pub alleles: Vec<Vec<u8>>,
    /// phred-scaled QUAL (-10 * log10PError)
    pub qual: f64,
    pub info: Vec<(&'static str, String)>,
    /// genotype allele indices; None = no-call (./.)
    pub gt: Option<(usize, usize)>,
    pub ad: Option<Vec<i32>>,
    pub dp: Option<i32>,
    pub gq: Option<i32>,
    pub pl: Vec<i32>,
    /// FORMAT SB (StrandBiasBySample, reference-confidence mode only)
    pub sb: Option<[i32; 4]>,
    /// physical phasing (reference-confidence mode only)
    pub phase: Option<Phase>,
    /// the other samples' genotypes (several samples: VCF mode, no SB / phasing)
    pub more: Vec<SampleGenotype>,
}

/// One more sample's genotype fields.
#[derive(Clone, Debug)]
pub struct SampleGenotype {
    pub gt: Option<(usize, usize)>,
    pub ad: Option<Vec<i32>>,
    pub dp: Option<i32>,
    pub gq: Option<i32>,
    pub pl: Vec<i32>,
}

/// `AssemblyBasedCallerUtils.phaseVC`: the genotype becomes phased, possibly with its alleles reversed,
/// and gets PGT / PID / PS.
#[derive(Clone, Debug)]
pub struct Phase {
    pub pgt: &'static str,
    pub pid: String,
    pub ps: i64,
    pub reversed: bool,
}

fn join_ints(v: &[i32]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

impl VcfRecord {
    /// The tab-separated data line (no trailing newline).
    pub fn to_line(&self) -> String {
        let mut s = String::with_capacity(256);
        s.push_str(&self.contig);
        s.push('\t');
        s.push_str(&self.pos.to_string());
        s.push('\t');
        s.push_str(self.id.as_deref().unwrap_or("."));
        s.push('\t');
        s.push_str(std::str::from_utf8(&self.alleles[0]).unwrap());
        s.push('\t');
        let alts: Vec<&str> = self.alleles[1..].iter().map(|a| std::str::from_utf8(a).unwrap()).collect();
        s.push_str(if alts.is_empty() { "." } else { "" });
        s.push_str(&alts.join(","));
        s.push('\t');
        s.push_str(&format_qual(self.qual));
        s.push_str("\t.\t");
        if self.info.is_empty() {
            s.push('.');
        } else {
            // VCFEncoder.writeInfoString: an empty value (a Flag) is written as the key alone
            let parts: Vec<String> = self.info.iter().map(|(k, v)| if v.is_empty() { k.to_string() } else { format!("{}={}", k, v) }).collect();
            s.push_str(&parts.join(";"));
        }
        if !self.more.is_empty() {
            self.push_genotypes(&mut s);
            return s;
        }
        // FORMAT: GT first, then the other keys present, sorted
        let mut keys = vec!["GT"];
        let sep = if self.phase.is_some() { '|' } else { '/' };
        let rev = self.phase.as_ref().map_or(false, |p| p.reversed);
        let mut vals = vec![match self.gt {
            Some((i, j)) if rev => format!("{}{}{}", j, sep, i),
            Some((i, j)) => format!("{}{}{}", i, sep, j),
            None => format!(".{}.", sep),
        }];
        if let Some(ad) = &self.ad {
            keys.push("AD");
            vals.push(join_ints(ad));
        }
        if let Some(dp) = self.dp {
            keys.push("DP");
            vals.push(dp.to_string());
        }
        if let Some(gq) = self.gq {
            keys.push("GQ");
            vals.push(gq.min(99).to_string()); // IntGenotypeFieldAccessors.GQAccessor caps at 99
        }
        if let Some(p) = &self.phase {
            keys.push("PGT");
            vals.push(p.pgt.to_string());
            keys.push("PID");
            vals.push(p.pid.clone());
        }
        if !self.pl.is_empty() {
            keys.push("PL");
            vals.push(join_ints(&self.pl));
        }
        if let Some(p) = &self.phase {
            keys.push("PS");
            vals.push(p.ps.to_string());
        }
        if let Some(sb) = &self.sb {
            keys.push("SB");
            vals.push(join_ints(sb));
        }
        s.push('\t');
        s.push_str(&keys.join(":"));
        s.push('\t');
        s.push_str(&vals.join(":"));
        s
    }

    /// FORMAT and sample columns of several samples (htsjdk `VCFEncoder`): GT, then the keys any genotype has,
    /// sorted; a missing value is ".", trailing missing values of a sample are dropped.
    fn push_genotypes(&self, s: &mut String) {
        let first = SampleGenotype { gt: self.gt, ad: self.ad.clone(), dp: self.dp, gq: self.gq, pl: self.pl.clone() };
        let all: Vec<&SampleGenotype> = std::iter::once(&first).chain(self.more.iter()).collect();
        let (has_ad, has_dp, has_gq, has_pl) = (
            all.iter().any(|g| g.ad.is_some()),
            all.iter().any(|g| g.dp.is_some()),
            all.iter().any(|g| g.gq.is_some()),
            all.iter().any(|g| !g.pl.is_empty()),
        );
        let mut keys = vec!["GT"];
        for (k, on) in [("AD", has_ad), ("DP", has_dp), ("GQ", has_gq), ("PL", has_pl)] {
            if on {
                keys.push(k);
            }
        }
        s.push('\t');
        s.push_str(&keys.join(":"));
        for g in all {
            let mut vals: Vec<String> = Vec::with_capacity(keys.len());
            if has_ad {
                vals.push(g.ad.as_ref().map_or(".".into(), |a| join_ints(a)));
            }
            if has_dp {
                vals.push(g.dp.map_or(".".into(), |d| d.to_string()));
            }
            if has_gq {
                vals.push(g.gq.map_or(".".into(), |q| q.min(99).to_string())); // GQAccessor caps at 99
            }
            if has_pl {
                vals.push(if g.pl.is_empty() { ".".into() } else { join_ints(&g.pl) });
            }
            while vals.last().is_some_and(|v| v == ".") {
                vals.pop();
            }
            s.push('\t');
            match g.gt {
                Some((i, j)) => s.push_str(&format!("{}/{}", i, j)),
                None => s.push_str("./."),
            }
            for v in vals {
                s.push(':');
                s.push_str(&v);
            }
        }
    }
}
