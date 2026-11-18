//! `--dbsnp`: the ID column and the `DB` flag of HaplotypeCaller calls (`VariantAnnotatorEngine.annotateContext` ->
//! `VariantOverlapAnnotator.annotateRsID` / `annotateOverlaps`). A dbSNP record starting at the call's position
//! (`FeatureContext.getValues(dbsnp, start)`), not filtered, matches when one of its biallelic events equals one of
//! the call's (`splitVariantContextToEvents` with left trimming, then `Event`'s minimal representation); the IDs of
//! the matching records, in file order, joined with ';', become the call's ID, and the call gets `DB`.

use crate::mutect::germline::Germline;
use hc_genotype::vcf::VcfRecord;
use std::collections::HashSet;

pub struct Dbsnp {
    db: Germline,
}

/// `Event`: start, ref, alt.
type Event = (i64, Vec<u8>, Vec<u8>);

/// `Allele.wouldBeSymbolicAllele`
fn is_symbolic(a: &[u8]) -> bool {
    a.len() > 1 && (a[0] == b'<' || a[a.len() - 1] == b'>' || a.contains(&b'[') || a.contains(&b']') || a[0] == b'.' || a[a.len() - 1] == b'.')
}

fn is_span_del(a: &[u8]) -> bool {
    a == b"*"
}

/// `GATKVariantContextUtils.trimAlleles(vc, trimForward = true, trimReverse = true)` on the alleles of one record.
fn trim_alleles(start: i64, alleles: &[Vec<u8>]) -> (i64, Vec<Vec<u8>>) {
    if alleles.len() <= 1 || alleles.iter().any(|a| a.len() == 1 && !is_span_del(a)) {
        return (start, alleles.to_vec());
    }
    let seqs: Vec<&Vec<u8>> = alleles.iter().filter(|a| !is_symbolic(a) && !is_span_del(a)).collect();
    // AlignmentUtils.normalizeAlleles(sequences, ranges, maxShift = 0, trim = true)
    let mut ranges: Vec<(usize, usize)> = seqs.iter().map(|s| (0, s.len())).collect();
    let mut min_size = ranges.iter().map(|r| r.1 - r.0).min().unwrap();
    let mut end_trim = 0usize;
    while min_size > 0 && seqs.iter().zip(&ranges).all(|(s, r)| s[r.1 - 1] == seqs[0][ranges[0].1 - 1]) {
        ranges.iter_mut().for_each(|r| r.1 -= 1);
        min_size -= 1;
        end_trim += 1;
    }
    let mut start_trim = 0usize;
    while min_size > 0 && seqs.iter().zip(&ranges).all(|(s, r)| s[r.0] == seqs[0][ranges[0].0]) {
        ranges.iter_mut().for_each(|r| r.0 += 1);
        min_size -= 1;
        start_trim += 1;
    }
    let empty = ranges.iter().any(|r| r.1 == r.0);
    let end_clip = if empty && start_trim == 0 { end_trim - 1 } else { end_trim };
    let start_clip = if empty && start_trim > 0 { start_trim - 1 } else { start_trim };
    if start_clip == 0 && end_clip == 0 {
        return (start, alleles.to_vec());
    }
    let out = alleles
        .iter()
        .map(|a| if is_symbolic(a) || is_span_del(a) { a.clone() } else { a[start_clip..a.len() - end_clip].to_vec() })
        .collect();
    (start + start_clip as i64, out)
}

/// `new Event(contig, start, ref, alt)`: identical suffixes removed (`makeMinimalRepresentation`).
fn event(start: i64, r: &[u8], a: &[u8]) -> Event {
    if r.len() == 1 || a.len() == 1 || r[r.len() - 1] != a[a.len() - 1] {
        return (start, r.to_vec(), a.to_vec());
    }
    assert!(r != a, "ref and alt alleles are identical");
    let min = r.len().min(a.len());
    let mut k = 0;
    while k < min && r[r.len() - 1 - k] == a[a.len() - 1 - k] {
        k += 1;
    }
    assert!(k < r.len() && k < a.len(), "Null alleles are not supported");
    (start, r[..r.len() - k].to_vec(), a[..a.len() - k].to_vec())
}

/// `splitVariantContextToEvents(vc, trimLeft = true, ...)`: biallelic records as they are, multiallelic ones split
/// per alt allele and trimmed.
fn events(start: i64, alleles: &[Vec<u8>]) -> Vec<Event> {
    match alleles.len() {
        0 | 1 => Vec::new(),
        2 => vec![event(start, &alleles[0], &alleles[1])],
        _ => alleles[1..]
            .iter()
            .map(|alt| {
                let (s, t) = trim_alleles(start, &[alleles[0].clone(), alt.clone()]);
                event(s, &t[0], &t[1])
            })
            .collect(),
    }
}

impl Dbsnp {
    pub fn load(path: &str, contigs: &HashSet<String>) -> Dbsnp {
        Dbsnp { db: Germline::load_with(path, contigs, false, "dbsnp") }
    }

    /// The IDs of the dbSNP records matching a call at `pos` with `alleles` (in file order).
    pub fn ids(&self, contig: &str, pos: i64, alleles: &[String]) -> Vec<String> {
        let al: Vec<Vec<u8>> = alleles.iter().map(|a| a.as_bytes().to_vec()).collect();
        self.db.with_records_at(contig, pos, |recs| {
            let mut ids: Vec<String> = Vec::new();
            if recs.is_empty() {
                return ids;
            }
            let mine = events(pos, &al);
            for r in recs {
                if r.filtered {
                    continue;
                }
                let alleles: Vec<Vec<u8>> = r.alleles.iter().map(|a| a.to_ascii_uppercase()).collect();
                if events(r.start, &alleles).iter().any(|e| mine.contains(e)) {
                    ids.push(r.id.clone());
                }
            }
            ids
        })
    }

    /// Sets the call's ID and `DB` flag from the dbSNP records starting at its position.
    pub fn annotate(&self, contig: &str, rec: &mut VcfRecord) {
        let ids = self.db.with_records_at(contig, rec.pos, |recs| {
            let mut ids: Vec<String> = Vec::new();
            if recs.is_empty() {
                return ids;
            }
            let mine = events(rec.pos, &rec.alleles);
            for r in recs {
                if r.filtered {
                    continue;
                }
                let alleles: Vec<Vec<u8>> = r.alleles.iter().map(|a| a.to_ascii_uppercase()).collect();
                if events(r.start, &alleles).iter().any(|e| mine.contains(e)) {
                    ids.push(r.id.clone());
                }
            }
            ids
        });
        if ids.is_empty() {
            return;
        }
        rec.id = Some(ids.join(";"));
        let at = rec.info.partition_point(|(k, _)| *k < "DB");
        rec.info.insert(at, ("DB", String::new()));
    }
}

/// `##INFO` line of the `DB` flag (`VCFStandardHeaderLines`).
pub const DB_HEADER: &str = "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">";

#[cfg(test)]
mod tests {
    use super::*;
    fn v(s: &[&str]) -> Vec<Vec<u8>> {
        s.iter().map(|x| x.as_bytes().to_vec()).collect()
    }
    #[test]
    fn events_like_gatk() {
        // biallelic: suffix trimmed by Event only
        assert_eq!(events(10, &v(&["ATT", "AT"])), vec![(10, b"AT".to_vec(), b"A".to_vec())]);
        assert_eq!(events(10, &v(&["AC", "AG"])), vec![(10, b"AC".to_vec(), b"AG".to_vec())]);
        // multiallelic: each pair trimmed on both sides
        assert_eq!(
            events(10, &v(&["AC", "AG", "A"])),
            vec![(11, b"C".to_vec(), b"G".to_vec()), (10, b"AC".to_vec(), b"A".to_vec())]
        );
        assert_eq!(
            events(10, &v(&["ATT", "AT", "A", "<NON_REF>"]))[..2],
            [(10, b"AT".to_vec(), b"A".to_vec()), (10, b"ATT".to_vec(), b"A".to_vec())]
        );
        // insertion pair padded on both sides: one base restored at the start
        assert_eq!(events(5, &v(&["GAT", "GAAT", "C"]))[0], (5, b"G".to_vec(), b"GA".to_vec()));
        assert!(events(5, &v(&["A"])).is_empty());
    }
}
