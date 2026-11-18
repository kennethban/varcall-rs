//! Output: VCF/GVCF headers and the final (bgzipped + tabix-indexed, or plain) files assembled from
//! the per-contig parts in output order.

use crate::backend::Part;
use crate::bgzf::{self, BgzfWriter};
use crate::tabix::{self, ContigIndex, VcfIndexer};
use std::fs::File;
use std::io::{self, BufWriter, Write};

pub const VCF_TEMPLATE: &str = include_str!("../tables/vcf_header_template.txt");
pub const GVCF_TEMPLATE: &str = include_str!("../tables/gvcf_header_template.txt");
pub const RAW_TEMPLATE: &str = include_str!("../tables/raw_header_template.txt");

/// `##INFO` lines `-G AS_StandardAnnotation` adds in GVCF mode (`getRawDescriptions` of the reducible
/// annotations, `getDescriptions` of AS_InbreedingCoeff and AS_QD), as GATK 4.7 writes them.
pub const AS_GVCF_HEADER: &[&str] = &[
    "##INFO=<ID=AS_InbreedingCoeff,Number=A,Type=Float,Description=\"Allele-specific inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">",
    "##INFO=<ID=AS_QD,Number=A,Type=Float,Description=\"Allele-specific Variant Confidence/Quality by Depth\">",
    "##INFO=<ID=AS_RAW_BaseQRankSum,Number=1,Type=String,Description=\"raw data for allele specific rank sum test of base qualities\">",
    "##INFO=<ID=AS_RAW_MQ,Number=1,Type=String,Description=\"Allele-specfic raw data for RMS Mapping Quality\">",
    "##INFO=<ID=AS_RAW_MQRankSum,Number=1,Type=String,Description=\"Allele-specfic raw data for Mapping Quality Rank Sum\">",
    "##INFO=<ID=AS_RAW_ReadPosRankSum,Number=1,Type=String,Description=\"allele specific raw data for rank sum test of read position bias\">",
    "##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.\">",
];

/// One `##contig` line (htsjdk `VCFContigHeaderLine(SAMSequenceRecord, assembly)`: ID, length, and
/// `assembly` when the BAM's `@SQ` line has `AS`).
pub struct ContigLine {
    pub name: String,
    pub len: i64,
    pub assembly: Option<String>,
}

fn header_value(v: &str) -> String {
    // VCFHeaderLine.toStringEncoding quotes values with commas or spaces
    if v.contains(',') || v.contains(' ') {
        format!("\"{}\"", v.replace('"', "\\\""))
    } else {
        v.to_string()
    }
}

/// The header template adjusted for this run's options: lines starting with one of `drop` removed, `extra` added,
/// metadata lines sorted as htsjdk's writer sorts them (`VCFHeader.getMetaDataInSortedOrder`: a `TreeSet` ordered by
/// each line's text); `##fileformat` stays first and `#CHROM` last.
pub fn template_with(template: &str, drop: &[&str], extra: &[String]) -> String {
    let mut first = None;
    let mut chrom = None;
    let mut meta: Vec<String> = Vec::new();
    for l in template.lines() {
        if l.starts_with("##fileformat=") {
            first = Some(l.to_string());
        } else if l.starts_with("#CHROM") {
            chrom = Some(l.to_string());
        } else if !drop.iter().any(|d| l.starts_with(d)) {
            meta.push(l.to_string());
        }
    }
    meta.extend(extra.iter().cloned());
    meta.sort();
    meta.dedup();
    let mut out = String::new();
    for l in first.iter().chain(meta.iter()).chain(chrom.iter()) {
        out.push_str(l);
        out.push('\n');
    }
    out
}

/// `GVCFBlockCombiner.rangeToVCFHeaderLine` for each band of `--gvcf-gq-bands`.
pub fn gvcf_block_lines(bands: &[i32]) -> Vec<String> {
    let mut out = Vec::new();
    let mut last = 0;
    for &v in bands.iter().chain(if *bands.last().unwrap() <= 99 { &[100][..] } else { &[][..] }) {
        out.push(format!("##GVCFBlock{last}-{v}=minGQ={last}(inclusive),maxGQ={v}(exclusive)"));
        last = v;
    }
    out
}

/// `GVCFBlockCombiner.parsePartitions` checks (GATK throws IllegalArgumentException).
pub fn check_gq_bands(bands: &[i32]) -> Result<(), String> {
    let mut last = 0;
    for &v in bands {
        if v < 0 {
            return Err("The list of GQ partitions contains a non-positive integer.".into());
        } else if v > 100 {
            return Err(format!("The value {v} in the list of GQ partitions is greater than VCFConstants.MAX_GENOTYPE_QUAL + 1 = 100."));
        } else if v < last {
            return Err(format!("The list of GQ partitions is out of order. Previous value is {last} but the next is {v}."));
        } else if v == last {
            return Err(format!("The value {v} appears more than once in the list of GQ partitions."));
        }
        last = v;
    }
    Ok(())
}

/// GATK's HaplotypeCaller header with this run's contigs and sample. Every line is GATK's, in GATK's
/// order, except `##GATKCommandLine`, which is replaced at the same position by
/// `##variant-rsCommandLine=<ID=variant-rs,CommandLine="...",Version="...",Date="...">`.
/// htsjdk's rule for a block-compressed (BGZF + .tbi) VCF output (IOUtil.hasBlockCompressedExtension): the name
/// ends, ignoring case, in .gz, .gzip, .bgz or .bgzf. Any other name is written as plain text.
pub fn block_compressed(path: &str) -> bool {
    let p = path.to_ascii_lowercase();
    [".gz", ".gzip", ".bgz", ".bgzf"].iter().any(|e| p.ends_with(e))
}

pub fn make_header(template: &str, contigs: &[ContigLine], samples: &[String], cmdline: &str, date: &str) -> String {
    let mut out = String::new();
    let mut contigs_done = false;
    let contig_lines = |out: &mut String| {
        for c in contigs {
            out.push_str(&format!("##contig=<ID={},length={}", header_value(&c.name), c.len));
            if let Some(a) = &c.assembly {
                out.push_str(&format!(",assembly={}", header_value(a)));
            }
            out.push_str(">\n");
        }
    };
    for l in template.lines() {
        if l.starts_with("##GATKCommandLine") {
            out.push_str(&format!(
                "##variant-rsCommandLine=<ID=variant-rs,CommandLine=\"{}\",Version=\"{}\",Date=\"{}\">\n",
                cmdline.replace('"', "\\\""),
                env!("CARGO_PKG_VERSION"),
                date
            ));
        } else if l.starts_with("##contig=") {
            if !contigs_done {
                contig_lines(&mut out);
                contigs_done = true;
            }
        } else if l.starts_with("#CHROM") {
            if !contigs_done {
                contig_lines(&mut out);
                contigs_done = true;
            }
            let mut f: Vec<&str> = l.split('\t').collect();
            f.pop();
            f.extend(samples.iter().map(|s| s.as_str()));
            out.push_str(&f.join("\t"));
            out.push('\n');
        } else {
            out.push_str(l);
            out.push('\n');
        }
    }
    out
}

/// A final output file.
pub struct OutFile {
    pub path: String,
    w: BufWriter<File>,
    offset: u64,
    level: Option<i32>,
    index: Vec<ContigIndex>,
    pub lines: u64,
    pub deflate_secs: f64,
}

impl OutFile {
    /// Creates the file and writes the header. `level` = Some(compression level) for bgzip + tabix.
    pub fn create(path: &str, header: &str, level: Option<i32>) -> io::Result<OutFile> {
        let w = BufWriter::with_capacity(4 << 20, File::create(path).map_err(|e| io::Error::new(e.kind(), format!("{path}: {e}")))?);
        let mut o = OutFile { path: path.to_string(), w, offset: 0, level, index: Vec::new(), lines: 0, deflate_secs: 0.0 };
        o.add_text(header.as_bytes())?;
        Ok(o)
    }

    /// Appends text (complete lines). Compressed into its own blocks.
    pub fn add_text(&mut self, text: &[u8]) -> io::Result<()> {
        if text.is_empty() {
            return Ok(());
        }
        match self.level {
            None => {
                self.w.write_all(text)?;
                self.offset += text.len() as u64;
            }
            Some(level) => {
                let mut bw = BgzfWriter::new(&mut self.w, level, VcfIndexer::default(), true);
                bw.write_all(text)?;
                let (_, sizes, mut idx, secs) = bw.finish()?;
                idx.finish();
                self.deflate_secs += secs;
                self.lines += idx.n_lines;
                let cum = bgzf::cumulative(self.offset, &sizes);
                for mut c in idx.contigs {
                    c.map_offsets(|v| bgzf::real_voff(v, &cum));
                    self.push_index(c);
                }
                self.offset = *cum.last().unwrap();
            }
        }
        Ok(())
    }

    /// Appends a finished part (same format as this file).
    pub fn add_part(&mut self, part: &Part) -> io::Result<()> {
        let mut f = File::open(&part.path)?;
        let n = io::copy(&mut f, &mut self.w)?;
        match (&part.block_sizes, self.level) {
            (Some(sizes), Some(_)) => {
                let cum = bgzf::cumulative(self.offset, sizes);
                assert_eq!(*cum.last().unwrap() - self.offset, n, "part size mismatch");
                for c in &part.index {
                    let mut c2 = c.clone();
                    c2.map_offsets(|v| bgzf::real_voff(v, &cum));
                    self.push_index(c2);
                }
            }
            (None, None) => {}
            _ => panic!("part format does not match the output format"),
        }
        self.offset += n;
        Ok(())
    }

    /// Adds a contig's index; consecutive pieces of one contig (split contigs) become one entry, as
    /// a .tbi must list each sequence once.
    fn push_index(&mut self, c: ContigIndex) {
        match self.index.last_mut() {
            Some(last) if last.name == c.name => last.merge(c),
            _ => self.index.push(c),
        }
    }

    /// Writes the EOF marker and the .tbi.
    pub fn finish(mut self) -> io::Result<()> {
        if self.level.is_some() {
            self.w.write_all(&bgzf::EOF_BLOCK)?;
        }
        self.w.flush()?;
        drop(self.w);
        if self.level.is_some() {
            tabix::write_tbi(&format!("{}.tbi", self.path), &self.index)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn templates_are_in_htsjdk_order() {
        for t in [VCF_TEMPLATE, GVCF_TEMPLATE] {
            let contigs = |x: &str| x.lines().filter(|l| !l.starts_with("##contig=")).map(|l| l.to_string()).collect::<Vec<_>>();
            assert_eq!(contigs(&template_with(t, &[], &[])), contigs(t));
        }
        let default: Vec<i32> = (1..=60).chain([70, 80, 90, 99]).collect();
        let lines: Vec<String> = GVCF_TEMPLATE.lines().filter(|l| l.starts_with("##GVCFBlock")).map(|l| l.to_string()).collect();
        let mut mine = gvcf_block_lines(&default);
        mine.sort();
        assert_eq!(mine, lines);
        assert!(check_gq_bands(&[10, 20, 30]).is_ok());
        assert!(check_gq_bands(&[0, 20]).is_err());
        assert!(check_gq_bands(&[20, 10]).is_err());
        assert!(check_gq_bands(&[101]).is_err());
    }
}
