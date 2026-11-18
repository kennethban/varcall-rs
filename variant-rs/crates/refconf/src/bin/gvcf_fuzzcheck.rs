//! Checks `GvcfWriter` (banded GVCF + unbanded records) against GATK's `GVCFWriter` on the synthetic record
//! streams written by `scratch/fuzz/GvcfFuzz.java` (a development probe, not in the repository).
//!
//! usage: gvcf_fuzzcheck <prefix>   (reads <prefix>.in, compares with <prefix>.gvcf and <prefix>.raw)

use hc_genotype::VcfRecord;
use hc_refconf::{default_gq_bands, CallRecord, GvcfWriter, Item, SiteRecord};

fn ints(s: &str) -> Vec<i32> {
    s.split(',').map(|x| x.parse().unwrap()).collect()
}

fn gt(s: &str) -> Option<(usize, usize)> {
    if s == "." {
        return None;
    }
    let (a, b) = s.split_once('/').unwrap();
    Some((a.parse().unwrap(), b.parse().unwrap()))
}

fn data_lines(p: &str) -> Vec<String> {
    std::fs::read_to_string(p).unwrap().lines().filter(|l| !l.starts_with('#')).map(|l| l.to_string()).collect()
}

fn main() {
    let prefix = std::env::args().nth(1).expect("usage: gvcf_fuzzcheck <prefix>");
    let input = std::fs::read_to_string(format!("{prefix}.in")).unwrap();
    let mut w: GvcfWriter<Vec<u8>, Vec<u8>> = GvcfWriter::new("chr1", Vec::new(), Some(Vec::new()), &default_gq_bands());
    let mut region: Vec<Item> = Vec::new();
    let mut started = false;
    for l in input.lines() {
        let t: Vec<&str> = l.split(' ').collect();
        match t[0] {
            "R" => {
                if started {
                    w.add_region(std::mem::take(&mut region)).unwrap();
                }
                started = true;
            }
            "S" => {
                let pl = ints(t[6]);
                region.push(Item::Site(SiteRecord {
                    pos: t[1].parse().unwrap(),
                    ref_base: t[2].as_bytes()[0],
                    ad: [t[3].parse().unwrap(), t[4].parse().unwrap()],
                    gq: t[5].parse().unwrap(),
                    pl: [pl[0], pl[1], pl[2]],
                }));
            }
            "C" => {
                // C pos end ref alt qual dp gt gq ad0 pl ad1 fdp 0
                let pos: i64 = t[1].parse().unwrap();
                let g = gt(t[7]);
                let called = g.is_some();
                let rec = VcfRecord {
                    contig: "chr1".into(),
                    pos,
                    id: None,
                    alleles: vec![t[3].as_bytes().to_vec(), t[4].as_bytes().to_vec(), b"<NON_REF>".to_vec()],
                    qual: t[5].parse().unwrap(),
                    info: vec![("DP", t[6].to_string())],
                    gt: g,
                    ad: if called { Some(vec![t[9].parse().unwrap(), t[11].parse().unwrap(), 0]) } else { None },
                    dp: if called { Some(t[12].parse().unwrap()) } else { None },
                    gq: Some(t[8].parse().unwrap()),
                    pl: ints(t[10]),
                    sb: None,
                    phase: None,
                    more: Vec::new(),
                };
                region.push(Item::Call(CallRecord { pos, end: t[2].parse().unwrap(), rec }));
            }
            "B" => {
                // B pos ref qual gt gq dp pl
                let pos: i64 = t[1].parse().unwrap();
                let dp: i32 = t[6].parse().unwrap();
                let rec = VcfRecord {
                    contig: "chr1".into(),
                    pos,
                    id: None,
                    alleles: vec![t[2].as_bytes().to_vec(), b"<NON_REF>".to_vec()],
                    qual: t[3].parse().unwrap(),
                    info: vec![("DP", "7".to_string())],
                    gt: gt(t[4]),
                    ad: None,
                    dp: if dp >= 0 { Some(dp) } else { None },
                    gq: Some(t[5].parse().unwrap()),
                    pl: ints(t[7]),
                    sb: None,
                    phase: None,
                    more: Vec::new(),
                };
                region.push(Item::Call(CallRecord { pos, end: pos, rec }));
            }
            x => panic!("bad record {x}"),
        }
    }
    if started {
        w.add_region(region).unwrap();
    }
    w.finish().unwrap();
    let mut ok = true;
    for (name, got, file) in [("banded", &w.out_buffer(), format!("{prefix}.gvcf")), ("raw", &w.raw_buffer(), format!("{prefix}.raw"))] {
        let exp = data_lines(&file);
        let got: Vec<String> = String::from_utf8_lossy(got).lines().map(|l| l.to_string()).collect();
        let mut diff = 0;
        for (i, (g, e)) in got.iter().zip(exp.iter()).enumerate() {
            if g != e {
                if diff < 5 {
                    println!("{name} line {i}:\n- GATK: {e}\n+ rust: {g}");
                }
                diff += 1;
            }
        }
        println!("{name}: GATK {} lines, rust {} lines, differing {}", exp.len(), got.len(), diff);
        ok &= diff == 0 && exp.len() == got.len();
    }
    println!("paths: {}", hc_refconf::stats::report());
    println!("RESULT: {}", if ok { "IDENTICAL" } else { "DIFFERENT" });
    if !ok {
        std::process::exit(1);
    }
}
