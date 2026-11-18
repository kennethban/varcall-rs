//! Validator for hc_downsample against an instrumented GATK HaplotypeCaller run.
//!
//! Inputs
//! * `--log`   DS_PROBE_LOG written by GATK run with the patched classes of `probes/` on the classpath;
//! * `--bam`   the BAM GATK read; the validator streams it with `samtools view -M` over the padded
//!             intervals and applies HaplotypeCaller's read filters itself (independent of GATK);
//! * `-L`      the (sorted, merged) intervals given to GATK, `contig:start-end`, repeatable;
//! * `--target` `--max-reads-per-alignment-start`;
//! * `--vcf`   (optional) GATK's VCF, to compare the predicted jittered QD text.
//!
//! Checks
//! 1. every draw from GATK's global Random comes from ReservoirDownsampler.submit or
//!    QualByDepth.fixTooHighQD;
//! 2. the filtered BAM stream equals the stream GATK submitted to its PositionalDownsampler;
//! 3. ORACLE: replaying GATK's own interleaving of group draws and QD draws, the Rust downsampler +
//!    JavaRandom give GATK's exact kept-read stream (set and order) and QD values;
//! 4. RULE: with the interleaving computed by the schedule rule from the BAM, the region pops and
//!    the per-region number of QD draws only, the same stream, QD values, RNG step count; plus the
//!    ReadyTracker's readiness loci equal GATK's.
use hc_downsample::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::process::{Command, Stdio};

#[derive(Clone, Debug)]
struct Rd {
    name: String,
    flag: u16,
    start: i64,
    span: Option<(i64, i64)>,
}
impl Rd {
    fn key(&self) -> (String, u16, i64) {
        (self.name.clone(), self.flag, self.start)
    }
}

#[derive(Default)]
struct Shard {
    intervals: Vec<(i64, i64)>,
    contig: String,
    ins: Vec<(String, u16, i64)>,
    outs: Vec<(String, u16, i64)>,
    /// tokens in log order: G(start) at first reservoir draw of a group, Q(value) per jittered QD
    tokens: Vec<Tok>,
    pops: Vec<(i64, i64, Option<i64>)>,
    readies: Vec<(i64, i64, Option<i64>)>,
    applies: Vec<(i64, i64, bool)>,
    /// QD values per apply index
    qd_per_apply: Vec<Vec<(i64, f64)>>,
    n_ds: u64,
    n_qd: u64,
}

#[derive(Clone, Debug)]
enum Tok {
    G(i64),
    Q(#[allow(dead_code)] i64, f64),
}

fn parse_ivl(s: &str) -> (String, i64, i64) {
    let (c, r) = s.rsplit_once(':').expect("contig:start-end");
    let (a, b) = r.split_once('-').expect("contig:start-end");
    (c.to_string(), a.replace(',', "").parse().unwrap(), b.replace(',', "").parse().unwrap())
}

fn op_len(cigar: &str) -> Vec<(u8, i64)> {
    let mut v = Vec::new();
    let mut n = 0i64;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + (b - b'0') as i64;
        } else {
            v.push((b, n));
            n = 0;
        }
    }
    v
}

/// HaplotypeCaller default read filters (same logic as variant-rs passes_filters).
fn passes(flag: u16, mapq: i64, rname: &str, pos: i64, cigar: &[(u8, i64)], seq: &str, qual: &str, has_rg: bool, ctg_len: i64) -> bool {
    if mapq < 20 || mapq == 255 { return false; }
    if flag & 0x4 != 0 || rname == "*" || pos < 1 { return false; }
    if flag & 0x100 != 0 || flag & 0x400 != 0 || flag & 0x200 != 0 { return false; }
    let reflen: i64 = cigar.iter().filter(|e| b"MDN=X".contains(&e.0)).map(|e| e.1).sum();
    if reflen == 0 { return false; }
    // CigarUtils.isGood
    let n = cigar.len();
    for (i, &(op, len)) in cigar.iter().enumerate() {
        if len == 0 { return false; }
        if op == b'H' && i != 0 && i != n - 1 { return false; }
        if op == b'S' {
            let l = cigar[..i].iter().all(|e| e.0 == b'H');
            let r = cigar[i + 1..].iter().all(|e| e.0 == b'H');
            if !l && !r { return false; }
        }
    }
    let mut prev = false;
    for &(op, _) in cigar {
        let ind = op == b'I' || op == b'D';
        if prev && ind { return false; }
        prev = ind;
    }
    for e in cigar.iter() { if e.0 == b'D' { return false; } if !(e.0 == b'S' || e.0 == b'H') { break; } }
    for e in cigar.iter().rev() { if e.0 == b'D' { return false; } if !(e.0 == b'S' || e.0 == b'H') { break; } }
    // WellformedReadFilter
    let readlen: i64 = cigar.iter().filter(|e| b"MIS=X".contains(&e.0)).map(|e| e.1).sum();
    has_rg && seq != "*" && qual != "*" && qual.len() == seq.len() && readlen == seq.len() as i64
        && !cigar.iter().any(|e| e.0 == b'N') && pos <= ctg_len
}

fn contig_lengths(samtools: &str, bam: &str) -> HashMap<String, i64> {
    let o = Command::new(samtools).args(["view", "-H", bam]).output().expect("samtools");
    let mut m = HashMap::new();
    for l in String::from_utf8_lossy(&o.stdout).lines() {
        if !l.starts_with("@SQ") { continue; }
        let (mut sn, mut ln) = (String::new(), 0i64);
        for f in l.split('\t') {
            if let Some(v) = f.strip_prefix("SN:") { sn = v.to_string(); }
            if let Some(v) = f.strip_prefix("LN:") { ln = v.parse().unwrap(); }
        }
        m.insert(sn, ln);
    }
    m
}

/// Filtered reads of one shard, in BAM order: reads overlapping the intervals padded by 100.
fn shard_reads(samtools: &str, bam: &str, contig: &str, ivls: &[(i64, i64)], ctg_len: i64) -> Vec<Rd> {
    let regions: Vec<String> = ivls.iter().map(|&(s, e)| format!("{}:{}-{}", contig, (s - 100).max(1), (e + 100).min(ctg_len))).collect();
    let mut child = Command::new(samtools).arg("view").arg("-M").arg(bam).args(&regions).stdout(Stdio::piped()).spawn().expect("samtools");
    let rd = BufReader::with_capacity(1 << 20, child.stdout.take().unwrap());
    let mut v = Vec::new();
    for line in rd.lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        let flag: u16 = f[1].parse().unwrap();
        let pos: i64 = f[3].parse().unwrap();
        let mapq: i64 = f[4].parse().unwrap();
        let cig = if f[5] == "*" { Vec::new() } else { op_len(f[5]) };
        let has_rg = f[11..].iter().any(|t| t.starts_with("RG:Z:"));
        if !passes(flag, mapq, f[2], pos, &cig, f[9], f[10], has_rg, ctg_len) { continue; }
        let reflen: i64 = cig.iter().filter(|e| b"MDN=X".contains(&e.0)).map(|e| e.1).sum();
        let end = pos + reflen - 1;
        let pnext: i64 = f[7].parse().unwrap();
        let tlen: i32 = f[8].parse().unwrap();
        let mate_start = if flag & 0x1 != 0 && flag & 0x8 == 0 && f[6] != "*" { pnext } else { 0 };
        v.push(Rd { name: f[0].to_string(), flag, start: pos, span: pileup_span(pos, end, flag, mate_start, tlen) });
    }
    child.wait().unwrap();
    v
}

fn parse_log(path: &str) -> (Vec<Shard>, HashMap<String, u64>) {
    let f = BufReader::with_capacity(1 << 20, std::fs::File::open(path).expect("log"));
    let mut shards: Vec<Shard> = Vec::new();
    let mut callers: HashMap<String, u64> = HashMap::new();
    let mut last_in_start = i64::MIN;
    let mut last_g = i64::MIN;
    for line in f.lines() {
        let line = line.unwrap();
        let mut it = line.splitn(2, ' ');
        let tag = it.next().unwrap();
        let rest = it.next().unwrap_or("");
        let w: Vec<&str> = rest.split(' ').collect();
        match tag {
            "NEWDS" => { shards.push(Shard::default()); last_in_start = i64::MIN; last_g = i64::MIN; }
            "SHARD" => {
                let s = shards.last_mut().unwrap();
                for iv in rest.trim_matches(|c| c == '[' || c == ']').split(", ") {
                    let (c, a, b) = parse_ivl(iv);
                    s.contig = c;
                    s.intervals.push((a, b));
                }
            }
            "IN" | "OUT" => {
                let k = (w[0].to_string(), w[1].parse().unwrap(), w[3].parse().unwrap());
                let s = shards.last_mut().unwrap();
                if tag == "IN" { last_in_start = k.2; s.ins.push(k); } else { s.outs.push(k); }
            }
            "N" => {
                let caller = w[2].to_string();
                *callers.entry(caller.clone()).or_default() += 1;
                let s = shards.last_mut().unwrap();
                if caller.ends_with("ReservoirDownsampler.submit") {
                    s.n_ds += 1;
                    if last_in_start != last_g {
                        s.tokens.push(Tok::G(last_in_start));
                        last_g = last_in_start;
                    }
                } else if caller.ends_with("QualByDepth.fixTooHighQD") {
                    s.n_qd += 1;
                }
            }
            "POP" | "READY" => {
                let l = if w[2] == "END" { None } else { Some(w[2].parse().unwrap()) };
                let t = (w[0].parse().unwrap(), w[1].parse().unwrap(), l);
                let s = shards.last_mut().unwrap();
                if tag == "POP" { s.pops.push(t) } else { s.readies.push(t) }
            }
            "APPLY" => {
                let s = shards.last_mut().unwrap();
                s.applies.push((w[1].parse().unwrap(), w[2].parse().unwrap(), w[3] == "true"));
                s.qd_per_apply.push(Vec::new());
            }
            "QDJ" => {
                let s = shards.last_mut().unwrap();
                let v: f64 = w[2].parse().unwrap();
                let pos: i64 = w[1].parse().unwrap();
                s.tokens.push(Tok::Q(pos, v));
                s.qd_per_apply.last_mut().expect("QDJ before APPLY").push((pos, v));
            }
            _ => {}
        }
    }
    (shards, callers)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut log = String::new();
    let mut bam = String::new();
    let mut ivls: Vec<(String, i64, i64)> = Vec::new();
    let mut target = 50usize;
    let mut vcf: Option<String> = None;
    let mut samtools = "samtools".to_string();
    let mut dump: Option<String> = None;
    // deliberately wrong rules, to show the checks discriminate: exact | C | P1 | P2 | none
    let mut variant = "exact".to_string();
    let mut i = 1;
    while i < a.len() {
        match a[i].as_str() {
            "--log" => { log = a[i + 1].clone(); i += 1; }
            "--bam" => { bam = a[i + 1].clone(); i += 1; }
            "-L" => { ivls.push(parse_ivl(&a[i + 1])); i += 1; }
            "--target" => { target = a[i + 1].parse().unwrap(); i += 1; }
            "--vcf" => { vcf = Some(a[i + 1].clone()); i += 1; }
            "--samtools" => { samtools = a[i + 1].clone(); i += 1; }
            "--dump-kept" => { dump = Some(a[i + 1].clone()); i += 1; }
            "--variant" => { variant = a[i + 1].clone(); i += 1; }
            x => panic!("unknown arg {x}"),
        }
        i += 1;
    }
    let (shards, callers) = parse_log(&log);
    let lens = contig_lengths(&samtools, &bam);
    let mut fails = 0u32;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("[{}] {}: {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok { fails += 1; }
    };

    // 1. RNG callers
    let other: Vec<_> = callers.iter().filter(|(c, _)| !c.ends_with("ReservoirDownsampler.submit") && !c.ends_with("QualByDepth.fixTooHighQD")).collect();
    let total_steps: u64 = callers.values().sum();
    check("rng callers", other.is_empty(), format!("{:?}", callers));

    // shards in log order; -L intervals grouped by contig must match
    let mut by_contig: Vec<(String, Vec<(i64, i64)>)> = Vec::new();
    for (c, s, e) in &ivls {
        match by_contig.last_mut() { Some((cc, v)) if cc == c => v.push((*s, *e)), _ => by_contig.push((c.clone(), vec![(*s, *e)])) }
    }
    check("shards", by_contig.len() == shards.len() && by_contig.iter().zip(&shards).all(|(b, s)| b.0 == s.contig && b.1 == s.intervals),
          format!("{} shard(s) in log, {} from -L", shards.len(), by_contig.len()));

    let mut rng_oracle = JavaRandom::gatk();
    let mut rng_rule = JavaRandom::gatk();
    let mut vcf_expect: Vec<(String, i64, f64)> = Vec::new();
    let mut dump_w = dump.map(|p| std::io::BufWriter::new(std::fs::File::create(p).unwrap()));
    for (si, sh) in shards.iter().enumerate() {
        let ctg_len = lens[&sh.contig];
        let reads = shard_reads(&samtools, &bam, &sh.contig, &sh.intervals, ctg_len);
        // 2. input stream
        let same_in = reads.len() == sh.ins.len() && reads.iter().zip(&sh.ins).all(|(r, k)| r.key() == *k);
        let ngroups = { let mut g = 0; let mut last = i64::MIN; for r in &reads { if r.start != last { g += 1; last = r.start; } } g };
        let mut sizes: HashMap<i64, usize> = HashMap::new();
        for r in &reads { *sizes.entry(r.start).or_default() += 1; }
        let ndraw_groups = sizes.values().filter(|&&n| n > target).count();
        let maxg = sizes.values().copied().max().unwrap_or(0);
        check(&format!("shard {} {}: filtered BAM stream == GATK downsampler input", si, sh.contig), same_in,
              format!("{} reads (GATK {}), {} start positions, {} with > {} reads (max {})", reads.len(), sh.ins.len(), ngroups, ndraw_groups, target, maxg));

        // groups from GATK's own input stream (isolates check 3 from the filter port)
        let mut groups: Vec<(i64, Vec<(String, u16, i64)>)> = Vec::new();
        for k in &sh.ins {
            match groups.last_mut() { Some((s, v)) if *s == k.2 => v.push(k.clone()), _ => groups.push((k.2, vec![k.clone()])) }
        }
        // 3. ORACLE interleaving
        let gidx: HashMap<i64, usize> = groups.iter().enumerate().map(|(i, g)| (g.0, i)).collect();
        let mut kept: Vec<Option<Vec<(String, u16, i64)>>> = groups.iter().map(|g| if g.1.len() <= target { Some(g.1.clone()) } else { None }).collect();
        let mut qd_ok = true;
        let mut nq = 0;
        for t in &sh.tokens {
            match t {
                Tok::G(p) => { let gi = gidx[p]; kept[gi] = Some(reservoir_downsample(groups[gi].1.clone(), target, &mut rng_oracle)); }
                Tok::Q(_, v) => { nq += 1; let x = rng_oracle.qd_jitter(); if x.to_bits() != v.to_bits() { qd_ok = false; } }
            }
        }
        let oracle_out: Vec<_> = kept.into_iter().flat_map(|k| k.expect("drawing group without logged draws")).collect();
        check(&format!("shard {} ORACLE kept stream == GATK", si), oracle_out == sh.outs, format!("{} kept of {} (GATK {})", oracle_out.len(), sh.ins.len(), sh.outs.len()));
        check(&format!("shard {} ORACLE QD values == GATK", si), qd_ok, format!("{} jittered QD", nq));

        // 4. RULE
        // 4a readiness from the pops
        let mut pops_at: HashMap<Option<i64>, Vec<i64>> = HashMap::new();
        for &(_, e, l) in &sh.pops { pops_at.entry(l).or_default().push((e + 100).min(ctg_len)); }
        let mut tr = ReadyTracker::new();
        let mut ready_at: Vec<Option<ReadyAt>> = vec![None; sh.pops.len()];
        for &(s, e) in &sh.intervals {
            for l in s..=e {
                if let Some(id) = tr.on_locus(l, pops_at.get(&Some(l)).map(|v| v.as_slice()).unwrap_or(&[])) { ready_at[id] = Some(ReadyAt::Locus(l)); }
            }
        }
        for id in tr.on_end(pops_at.get(&None).map(|v| v.as_slice()).unwrap_or(&[])) { ready_at[id] = Some(ReadyAt::End); }
        let gatk_ready: Vec<Option<ReadyAt>> = sh.readies.iter().map(|r| Some(match r.2 { Some(l) => ReadyAt::Locus(l), None => ReadyAt::End })).collect();
        let order_ok = sh.pops.iter().map(|p| (p.0, p.1)).eq(sh.applies.iter().map(|a| (a.0, a.1)));
        check(&format!("shard {} ReadyTracker == GATK READY loci", si), ready_at == gatk_ready && order_ok,
              format!("{} regions, {} ready at END", ready_at.len(), ready_at.iter().filter(|r| **r == Some(ReadyAt::End)).count()));

        // 4b sequential schedule
        let mut ds: PositionalDownsampler<Rd> = PositionalDownsampler::new(target);
        for r in &reads { ds.push(r.start, r.clone()); }
        ds.finish();
        let mut cov = LibsCoverage::new();
        let mut out: Vec<Rd> = Vec::new();
        let mut buf: Vec<Rd> = Vec::new();
        let mut rule_qd_ok = true;
        let mut n_qd_rule = 0;
        let mut horizon_hist: HashMap<&'static str, u64> = HashMap::new();
        let nreg = sh.applies.len();
        for i in 0..nreg {
            let next_ready = if i + 1 < nreg { ready_at[i + 1] } else { None };
            let n = lookahead_locus(next_ready, &sh.intervals);
            let c = match n {
                None => None,
                Some(n) => loop {
                    let lim = ds.first_pending_start().map(|s| s - 1).unwrap_or(i64::MAX);
                    if let Some(c) = cov.first_covered(n, lim) { break Some(c); }
                    if ds.first_pending_start().is_none() { break None; }
                    buf.clear();
                    ds.finalize_head(&mut rng_rule, &mut buf);
                    for r in buf.drain(..) { cov.add(r.span); out.push(r); }
                },
            };
            let h = match (variant.as_str(), c) {
                ("exact", _) | (_, None) => ds.horizon(c),
                ("C", Some(c)) => Horizon::Through(c),
                ("P1", Some(c)) => match ds.horizon(Some(c)) { Horizon::Through(_) => { let p = reads.iter().map(|r| r.start).find(|&s| s > c); p.map(Horizon::Through).unwrap_or(Horizon::All) } x => x },
                ("P2", Some(c)) => { let mut it = reads.iter().map(|r| r.start).filter(|&s| s > c); let p1 = it.next(); let p2 = p1.and_then(|p1| it.find(|&s| s > p1)); p2.map(Horizon::Through).unwrap_or(Horizon::All) }
                ("none", _) => Horizon::All,
                _ => panic!("unknown variant"),
            };
            *horizon_hist.entry(match (n, c, h) {
                (None, _, _) => "stream drained (R_i+1 absent/END/last locus)",
                (Some(_), None, _) => "no pileup after N",
                (Some(n), Some(c), _) if c > n => "C > N",
                _ => "C == N",
            }).or_default() += 1;
            buf.clear();
            match h {
                Horizon::Through(p) => ds.finalize_upto(p, &mut rng_rule, &mut buf),
                Horizon::All => ds.finalize_all(&mut rng_rule, &mut buf),
                Horizon::NeedInput => unreachable!(),
            }
            for r in buf.drain(..) { cov.add(r.span); out.push(r); }
            for &(pos, v) in &sh.qd_per_apply[i] {
                let x = rng_rule.qd_jitter();
                n_qd_rule += 1;
                if x.to_bits() != v.to_bits() { rule_qd_ok = false; if n_qd_rule < 20 { eprintln!("  QD mismatch at {}:{} predicted {} GATK {}", sh.contig, pos, x, v); } }
                vcf_expect.push((sh.contig.clone(), pos, x));
            }
        }
        buf.clear();
        ds.finalize_all(&mut rng_rule, &mut buf);
        out.extend(buf.drain(..));
        let rule_keys: Vec<_> = out.iter().map(|r| r.key()).collect();
        if let Some(w) = dump_w.as_mut() {
            use std::io::Write;
            for r in &out { writeln!(w, "{}\t{}\t{}\t{}", sh.contig, r.start, r.name, r.flag).unwrap(); }
        }
        let first_diff = rule_keys.iter().zip(&sh.outs).position(|(a, b)| a != b);
        check(&format!("shard {} RULE kept stream == GATK", si), rule_keys == sh.outs,
              format!("{} kept (GATK {}){}", rule_keys.len(), sh.outs.len(), first_diff.map(|d| format!(", first diff at #{d}: {:?} vs {:?}", rule_keys[d], sh.outs[d])).unwrap_or_default()));
        check(&format!("shard {} RULE QD values == GATK", si), rule_qd_ok, format!("{} jittered QD, {} drawing groups, {} nextInt draws; horizons: {:?}", n_qd_rule, ds.stats.drawing_groups, ds.stats.draws, horizon_hist));
    }
    check("RNG LCG steps (oracle / rule / GATK)", rng_oracle.steps == total_steps && rng_rule.steps == total_steps,
          format!("{} / {} / {}", rng_oracle.steps, rng_rule.steps, total_steps));

    // QD text in the VCF
    if let Some(v) = vcf {
        // raw text (bcftools would re-format the floats)
        let o = Command::new("bgzip").args(["-dc", &v]).output().expect("bgzip");
        let mut m: HashMap<(String, i64), Vec<String>> = HashMap::new();
        for l in String::from_utf8_lossy(&o.stdout).lines() {
            if l.starts_with('#') { continue; }
            let f: Vec<&str> = l.split('\t').collect();
            if let Some(q) = f[7].split(';').find_map(|kv| kv.strip_prefix("QD=")) {
                m.entry((f[0].to_string(), f[1].parse().unwrap())).or_default().push(q.to_string());
            }
        }
        let mut bad = 0;
        for (c, p, x) in &vcf_expect {
            let s = javarandom::format_fixed(*x, 2);
            if !m.get(&(c.clone(), *p)).map(|v| v.contains(&s)).unwrap_or(false) {
                bad += 1;
                if bad < 10 { eprintln!("  VCF QD mismatch {}:{} predicted {} have {:?}", c, p, s, m.get(&(c.clone(), *p))); }
            }
        }
        check("VCF QD text", bad == 0, format!("{} predicted jittered QD values found in the VCF, {} missing", vcf_expect.len() - bad, bad));
    }
    println!("{}", if fails == 0 { "ALL PASS" } else { "FAILURES" });
    std::process::exit(if fails == 0 { 0 } else { 1 });
}
