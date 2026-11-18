//! Sequence dictionary, `-L` parsing (strings and interval files) and GATK's interval handling:
//! sorting in dictionary order, merging of overlapping and adjacent intervals
//! (`IntervalMergingRule.ALL`), grouping by contig and padding (`IntervalUtils.getIntervalsWithFlanks`).

use std::collections::HashMap;

#[derive(Clone)]
pub struct Dict {
    pub names: Vec<String>,
    pub lens: Vec<i64>,
    pub index: HashMap<String, usize>,
}

impl Dict {
    pub fn new(names: Vec<String>, lens: Vec<i64>) -> Dict {
        let index = names.iter().enumerate().map(|(i, n)| (n.clone(), i)).collect();
        Dict { names, lens, index }
    }

    /// The reference dictionary: `<fasta minus extension>.dict` if present, else the `.fai`.
    pub fn from_reference(fasta: &str) -> Dict {
        let stem = fasta.strip_suffix(".gz").unwrap_or(fasta);
        let stem = stem.rsplit_once('.').map(|(a, _)| a).unwrap_or(stem);
        for d in [format!("{stem}.dict"), format!("{fasta}.dict")] {
            if let Ok(t) = std::fs::read_to_string(&d) {
                let (mut names, mut lens) = (vec![], vec![]);
                for l in t.lines().filter(|l| l.starts_with("@SQ")) {
                    let (mut sn, mut ln) = (None, None);
                    for f in l.split('\t') {
                        if let Some(v) = f.strip_prefix("SN:") {
                            sn = Some(v.to_string());
                        }
                        if let Some(v) = f.strip_prefix("LN:") {
                            ln = Some(v.parse::<i64>().unwrap());
                        }
                    }
                    names.push(sn.expect("SN in .dict"));
                    lens.push(ln.expect("LN in .dict"));
                }
                return Dict::new(names, lens);
            }
        }
        let fai = std::fs::read_to_string(format!("{fasta}.fai")).expect("reference .dict or .fai");
        let (mut names, mut lens) = (vec![], vec![]);
        for l in fai.lines() {
            let f: Vec<&str> = l.split('\t').collect();
            names.push(f[0].to_string());
            lens.push(f[1].parse().unwrap());
        }
        Dict::new(names, lens)
    }
}

/// One interval, 1-based inclusive, on contig `tid` of the dictionary.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Interval {
    pub tid: usize,
    pub start: i64,
    pub end: i64,
}

fn parse_pos(s: &str) -> Result<i64, String> {
    s.replace(',', "").parse::<i64>().map_err(|_| format!("bad position '{s}'"))
}

/// GATK interval string: `contig`, `contig:start`, `contig:start-end`, `contig:start+`. Contig names may
/// contain ':' (HLA alleles); the dictionary decides.
pub fn parse_interval(s: &str, dict: &Dict) -> Result<Interval, String> {
    let s = s.trim();
    let whole = dict.index.get(s).map(|&t| Interval { tid: t, start: 1, end: dict.lens[t] });
    let with_pos = s.rfind(':').and_then(|i| {
        let (c, p) = (&s[..i], &s[i + 1..]);
        let &t = dict.index.get(c)?;
        let len = dict.lens[t];
        let r = if let Some(a) = p.strip_suffix('+') {
            parse_pos(a).map(|a| (a, len))
        } else if let Some((a, b)) = p.split_once('-') {
            parse_pos(a).and_then(|a| parse_pos(b).map(|b| (a, b)))
        } else {
            parse_pos(p).map(|a| (a, a))
        };
        Some(r.map(|(a, b)| Interval { tid: t, start: a, end: b }))
    });
    let iv = match (whole, with_pos) {
        (Some(_), Some(Ok(_))) => return Err(format!("ambiguous interval '{s}': both a contig name and contig:position")),
        (Some(w), _) => w,
        (None, Some(r)) => r?,
        (None, None) => return Err(format!("interval '{s}': contig not in the sequence dictionary")),
    };
    let len = dict.lens[iv.tid];
    if iv.start < 1 || iv.end < iv.start || iv.end > len {
        return Err(format!("interval '{s}' is outside contig {} (length {len}) or empty", dict.names[iv.tid]));
    }
    Ok(iv)
}

/// Reads an interval file: `.bed` (0-based half-open), `.interval_list` (Picard: `@` header, then
/// contig, start, end, strand, name), anything else one GATK interval string per line.
pub fn read_interval_file(path: &str, dict: &Dict) -> Result<Vec<Interval>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut out = Vec::new();
    for l in text.lines() {
        let l = l.trim();
        if l.is_empty() || l.starts_with('#') || l.starts_with('@') || l.starts_with("track") || l.starts_with("browser") {
            continue;
        }
        if path.ends_with(".bed") {
            let f: Vec<&str> = l.split('\t').collect();
            let &t = dict.index.get(f[0]).ok_or(format!("{path}: unknown contig {}", f[0]))?;
            let (s, e) = (parse_pos(f[1])? + 1, parse_pos(f[2])?);
            if e >= s {
                out.push(Interval { tid: t, start: s, end: e });
            }
        } else if path.ends_with(".interval_list") {
            let f: Vec<&str> = l.split('\t').collect();
            let &t = dict.index.get(f[0]).ok_or(format!("{path}: unknown contig {}", f[0]))?;
            out.push(Interval { tid: t, start: parse_pos(f[1])?, end: parse_pos(f[2])? });
        } else {
            out.push(parse_interval(l, dict)?);
        }
    }
    Ok(out)
}

/// `-L` value: an existing file, or an interval string.
pub fn parse_l_arg(arg: &str, dict: &Dict) -> Result<Vec<Interval>, String> {
    let is_file = std::path::Path::new(arg).is_file()
        && [".intervals", ".interval_list", ".list", ".bed"].iter().any(|e| arg.ends_with(e));
    if is_file {
        read_interval_file(arg, dict)
    } else {
        parse_interval(arg, dict).map(|i| vec![i])
    }
}

/// Sorts in dictionary order and merges overlapping or abutting intervals.
pub fn sort_merge(mut v: Vec<Interval>) -> Vec<Interval> {
    v.sort();
    let mut out: Vec<Interval> = Vec::with_capacity(v.len());
    for i in v {
        match out.last_mut() {
            Some(l) if l.tid == i.tid && i.start <= l.end + 1 => l.end = l.end.max(i.end),
            _ => out.push(i),
        }
    }
    out
}

/// `IntervalUtils.loadIntervals` padding (`--interval-padding`, `--interval-exclusion-padding`): every interval
/// grows by `pad` on both sides, clipped to its contig (`createPaddedGenomeLoc`); merge afterwards.
pub fn pad(v: &mut [Interval], pad: i64, dict: &Dict) {
    if pad > 0 {
        for i in v {
            i.start = (i.start - pad).max(1);
            i.end = (i.end + pad).min(dict.lens[i.tid]);
        }
    }
}

/// `GenomeLocSortedSet.subtractRegions`: `include` minus `exclude` (both sorted and merged).
pub fn subtract(include: &[Interval], exclude: &[Interval]) -> Vec<Interval> {
    let mut out = Vec::new();
    let mut k = 0;
    for &i in include {
        let mut cur = i;
        let mut alive = true;
        while k < exclude.len() && (exclude[k].tid, exclude[k].end) < (i.tid, i.start) {
            k += 1;
        }
        let mut j = k;
        while alive && j < exclude.len() && exclude[j].tid == cur.tid && exclude[j].start <= cur.end {
            let x = exclude[j];
            if x.start > cur.start {
                out.push(Interval { tid: cur.tid, start: cur.start, end: x.start - 1 });
            }
            if x.end >= cur.end {
                alive = false;
            } else {
                cur.start = cur.start.max(x.end + 1);
            }
            j += 1;
        }
        if alive {
            out.push(cur);
        }
    }
    out
}

/// `IntervalUtils.getIntervalsWithFlanks(intervals, padding)`: pad, clip to the contig, merge.
pub fn pad_merge(v: &[(i64, i64)], pad: i64, len: i64) -> Vec<(i64, i64)> {
    let mut out: Vec<(i64, i64)> = Vec::new();
    for &(s, e) in v {
        let (s, e) = ((s - pad).max(1), (e + pad).min(len));
        match out.last_mut() {
            Some(l) if s <= l.1 + 1 => l.1 = l.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn subtract_regions() {
        let iv = |t, s, e| Interval { tid: t, start: s, end: e };
        let inc = vec![iv(0, 1, 100), iv(0, 200, 300), iv(1, 1, 50)];
        let exc = vec![iv(0, 10, 20), iv(0, 90, 210), iv(0, 250, 250), iv(1, 1, 50)];
        assert_eq!(subtract(&inc, &exc), vec![iv(0, 1, 9), iv(0, 21, 89), iv(0, 211, 249), iv(0, 251, 300)]);
        assert_eq!(subtract(&inc, &[]), inc);
    }
}

/// Groups sorted, merged intervals by contig (GATK `makeReadShards`: one read shard per contig).
pub fn group_by_contig(v: &[Interval]) -> Vec<(usize, Vec<(i64, i64)>)> {
    let mut out: Vec<(usize, Vec<(i64, i64)>)> = Vec::new();
    for i in v {
        match out.last_mut() {
            Some((t, l)) if *t == i.tid => l.push((i.start, i.end)),
            _ => out.push((i.tid, vec![(i.start, i.end)])),
        }
    }
    out
}
