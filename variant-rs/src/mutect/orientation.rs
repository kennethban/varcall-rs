//! Read-orientation (F1R2) artifact model: GATK 4.7 `F1R2CountsCollector` (Mutect2 `--f1r2-tar-gz`, fed from
//! the isActive tumor pileup), `LearnReadOrientationModel` (`variant-rs learnreadorientationmodel`) and the
//! artifact-prior tables read by FilterMutectCalls' ReadOrientationFilter (`--ob-priors`).
//!
//! Single tumor sample. Tables / metrics files are written as GATK writes them (the alt depth-1 histogram
//! columns are in a fixed order; GATK's order there depends on enum identity hashes).

use super::fmc_math::{beta_binomial_log_probability, java_double_to_string, jlog};
use super::somatic::normalize_from_log_to_linear;
use std::io::{Read, Write};

pub const MAX_DEPTH: usize = 200;
const MIN_MEDIAN_MQ: f64 = 50.0;

// ------------------------------------------------------------------------------------------------
// k-mers, states

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// `F1R2FilterConstants.ALL_KMERS`: htsjdk `generateAllKmers(3)` (AAC .. TTT, then AAA).
pub fn all_kmers() -> Vec<String> {
    let mut v: Vec<String> = Vec::with_capacity(64);
    for a in BASES {
        for b in BASES {
            for c in BASES {
                v.push(String::from_utf8(vec![a, b, c]).unwrap());
            }
        }
    }
    v.rotate_left(1);
    v
}

/// Lexicographic index of an ACGT k-mer (None if it has another character).
fn kmer_index(k: &[u8]) -> Option<usize> {
    let mut i = 0;
    for &b in k {
        i = i * 4 + base_index(b)?;
    }
    Some(i)
}

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn kmer_of(i: usize) -> String {
    String::from_utf8(vec![BASES[i / 16], BASES[(i / 4) % 4], BASES[i % 4]]).unwrap()
}

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            x => x as char,
        })
        .collect()
}

/// `F1R2FilterConstants.CANONICAL_KMERS` (first-appearance order over ALL_KMERS)
fn canonical_kmers() -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    for k in all_kmers() {
        let rc = revcomp(&k);
        let c = if k < rc { k } else { rc };
        if !v.contains(&c) {
            v.push(c);
        }
    }
    v
}

/// The ALL_KMERS-filled `HashMap<String, _>` iteration order (ref histograms, artifact prior tables).
fn kmer_hashmap_order() -> Vec<String> {
    let k = all_kmers();
    let refs: Vec<&str> = k.iter().map(|s| s.as_str()).collect();
    super::somatic::java_hashmap_order(&refs).into_iter().map(|i| k[i].clone()).collect()
}

/// `IOUtils.urlEncode` (java.net.URLEncoder, UTF-8)
pub fn url_encode(s: &str) -> String {
    let mut o = String::new();
    for b in s.bytes() {
        match b {
            b'a'..=b'z' | b'A'..=b'Z' | b'0'..=b'9' | b'.' | b'-' | b'*' | b'_' => o.push(b as char),
            b' ' => o.push('+'),
            _ => o.push_str(&format!("%{b:02X}")),
        }
    }
    o
}

// ArtifactState ordinals: F1R2_A..T (0-3), F2R1_A..T (4-7), HOM_REF 8, GERMLINE_HET 9, SOMATIC_HET 10, HOM_VAR 11
pub const NUM_STATES: usize = 12;
const HOM_REF: usize = 8;

fn rev_comp_state(s: usize) -> usize {
    match s {
        0..=3 => 4 + (3 - s), // F1R2_X -> F2R1_comp(X)
        4..=7 => 3 - (s - 4),
        _ => s,
    }
}

/// (alpha, beta) of the allele-fraction and alt-F1R2-fraction beta-binomials per state
fn pseudocounts(s: usize) -> ((f64, f64), (f64, f64)) {
    let af = match s {
        0..=7 => (1.0, 9.0),
        8 => (3.0, 10000.0),
        9 => (5.0, 5.0),
        10 => (2.0, 5.0),
        _ => (10000.0, 3.0),
    };
    let f1r2 = match s {
        0..=3 => (100.0, 1.0),
        4..=7 => (1.0, 100.0),
        _ => (10.0, 10.0),
    };
    (af, f1r2)
}

/// `LearnReadOrientationModelEngine.getFlatPrior`
fn flat_prior(ref_allele: usize) -> [f64; NUM_STATES] {
    let mut p = [1.0 / (NUM_STATES - 2) as f64; NUM_STATES];
    p[ref_allele] = 0.0;
    p[4 + ref_allele] = 0.0;
    p
}

/// `LearnReadOrientationModelEngine.computeResponsibilities` (bases as 0..3)
pub fn compute_responsibilities(ref_allele: usize, alt_allele: usize, alt_depth: i32, f1r2_alt: i32, depth: i32, prior: &[f64], given_not_hom_ref: bool) -> [f64; NUM_STATES] {
    let mut l = [0.0f64; NUM_STATES];
    for s in 0..NUM_STATES {
        if s == ref_allele || s == 4 + ref_allele {
            l[s] = f64::NEG_INFINITY;
            continue;
        }
        if s < 8 && s % 4 != alt_allele {
            l[s] = f64::NEG_INFINITY;
            continue;
        }
        let ((aa, ab), (fa, fb)) = pseudocounts(s);
        l[s] = jlog(prior[s]) + beta_binomial_log_probability(aa, ab, depth, alt_depth) + beta_binomial_log_probability(fa, fb, alt_depth, f1r2_alt);
    }
    if given_not_hom_ref {
        l[HOM_REF] = f64::NEG_INFINITY;
    }
    normalize_from_log_to_linear(&mut l);
    l
}

// ------------------------------------------------------------------------------------------------
// tar (ustar write; read with GNU long names / pax paths)

fn tar_header(name: &str, size: usize) -> [u8; 512] {
    let mut h = [0u8; 512];
    let nb = name.as_bytes();
    assert!(nb.len() < 100, "tar entry name too long: {name}");
    h[..nb.len()].copy_from_slice(nb);
    let oct = |h: &mut [u8; 512], off: usize, len: usize, v: u64| {
        let s = format!("{:0w$o}\0", v, w = len - 1);
        h[off..off + len].copy_from_slice(s.as_bytes());
    };
    oct(&mut h, 100, 8, 0o644);
    oct(&mut h, 108, 8, 0);
    oct(&mut h, 116, 8, 0);
    oct(&mut h, 124, 12, size as u64);
    let mtime = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs());
    oct(&mut h, 136, 12, mtime);
    h[156] = b'0';
    h[257..263].copy_from_slice(b"ustar\0");
    h[263..265].copy_from_slice(b"00");
    for b in &mut h[148..156] {
        *b = b' ';
    }
    let sum: u32 = h.iter().map(|&b| b as u32).sum();
    let s = format!("{:06o}\0 ", sum);
    h[148..156].copy_from_slice(s.as_bytes());
    h
}

pub fn write_tar_gz(path: &str, files: &[(String, Vec<u8>)]) {
    let f = std::fs::File::create(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut gz = flate2::write::GzEncoder::new(std::io::BufWriter::new(f), flate2::Compression::default());
    for (name, data) in files {
        // GATK's archives name their entries "./<file>" (tar of the output directory)
        gz.write_all(&tar_header(&format!("./{name}"), data.len())).unwrap();
        gz.write_all(data).unwrap();
        let pad = (512 - data.len() % 512) % 512;
        gz.write_all(&vec![0u8; pad]).unwrap();
    }
    gz.write_all(&[0u8; 1024]).unwrap();
    gz.finish().unwrap().flush().unwrap();
}

/// Regular files of a .tar.gz as (base name, contents), in archive order.
pub fn read_tar_gz(path: &str) -> Vec<(String, Vec<u8>)> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut data = Vec::new();
    flate2::read::MultiGzDecoder::new(f).read_to_end(&mut data).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut out = Vec::new();
    let mut off = 0;
    let mut long_name: Option<String> = None;
    let cstr = |b: &[u8]| String::from_utf8_lossy(&b[..b.iter().position(|&c| c == 0).unwrap_or(b.len())]).to_string();
    while off + 512 <= data.len() {
        let h = &data[off..off + 512];
        if h.iter().all(|&b| b == 0) {
            break;
        }
        let size_s = cstr(&h[124..136]);
        let size = usize::from_str_radix(size_s.trim(), 8).unwrap_or(0);
        let typ = h[156];
        let mut name = cstr(&h[0..100]);
        if &h[257..262] == b"ustar" {
            let prefix = cstr(&h[345..500]);
            if !prefix.is_empty() {
                name = format!("{prefix}/{name}");
            }
        }
        let body = &data[off + 512..off + 512 + size];
        off += 512 + size.div_ceil(512) * 512;
        match typ {
            b'L' => long_name = Some(cstr(body)),
            b'x' => {
                for rec in String::from_utf8_lossy(body).lines() {
                    if let Some(p) = rec.split_once(' ').and_then(|(_, kv)| kv.strip_prefix("path=")) {
                        long_name = Some(p.to_string());
                    }
                }
            }
            b'0' | 0 => {
                let n = long_name.take().unwrap_or(name);
                let base = n.rsplit('/').next().unwrap().to_string();
                out.push((base, body.to_vec()));
            }
            _ => {
                long_name = None;
            }
        }
    }
    out
}

// ------------------------------------------------------------------------------------------------
// F1R2CountsCollector

/// What F1R2CountsCollector needs from one locus' tumor pileup, restricted to elements with qual > 20 (deletions
/// never qualify: their qual is 16). Kept as counters, updated as reads are scattered and as overlap cleaning lowers
/// qualities before the locus is evaluated.
#[derive(Clone, Default)]
pub struct LocusCounts {
    /// elements with qual > 20 (including non-ACGT bases)
    pub n_kept: u32,
    pub cnt: [u32; 4],
    /// F1R2 reads per base
    pub f1r2: [u32; 4],
    /// isAfterInsertion || isBeforeDeletionStart
    pub indel_adj: u32,
    /// elements with mapping quality 60 (the usual case), and a (value, count) histogram of the others
    pub mq60: u32,
    pub mq: Vec<(u8, u32)>,
}

impl LocusCounts {
    pub fn clear(&mut self) {
        self.n_kept = 0;
        self.cnt = [0; 4];
        self.f1r2 = [0; 4];
        self.indel_adj = 0;
        self.mq60 = 0;
        self.mq.clear();
    }

    /// Adds (or removes) one element.
    #[inline]
    pub fn update(&mut self, base: u8, mapq: u8, f1r2: bool, indel_adj: bool, add: bool) {
        let d = |x: &mut u32| if add { *x += 1 } else { *x -= 1 };
        d(&mut self.n_kept);
        if let Some(b) = base_index(base) {
            d(&mut self.cnt[b]);
            if f1r2 {
                d(&mut self.f1r2[b]);
            }
        }
        if indel_adj {
            d(&mut self.indel_adj);
        }
        if mapq == 60 {
            d(&mut self.mq60);
            return;
        }
        match self.mq.iter().position(|e| e.0 == mapq) {
            Some(k) => {
                d(&mut self.mq[k].1);
                if self.mq[k].1 == 0 {
                    self.mq.swap_remove(k);
                }
            }
            None => {
                assert!(add, "removing an element that was never added");
                self.mq.push((mapq, 1));
            }
        }
    }

    /// `MathUtils.median(mappingQualities)` (commons-math Percentile 50, LEGACY) of the kept elements.
    fn median_mq(&self) -> f64 {
        if self.mq.is_empty() {
            return 60.0;
        }
        let mut h = self.mq.clone();
        if self.mq60 > 0 {
            h.push((60, self.mq60));
        }
        h.sort();
        let n = self.n_kept as usize;
        let at = |i: usize| -> f64 {
            let mut c = 0usize;
            for &(v, k) in &h {
                c += k as usize;
                if i < c {
                    return v as f64;
                }
            }
            unreachable!()
        };
        if n == 1 {
            return at(0);
        }
        let pos = 0.5 * (n as f64 + 1.0);
        let fpos = pos.floor();
        let ip = fpos as usize;
        let dif = pos - fpos;
        if pos >= n as f64 {
            return at(n - 1);
        }
        let (lower, upper) = (at(ip - 1), at(ip));
        lower + dif * (upper - lower)
    }
}

pub struct F1R2Counts {
    pub sample: String,
    /// [kmer index][depth - 1]
    ref_hist: Vec<[u64; MAX_DEPTH]>,
    /// [kmer index][alt base][orientation: 0 F1R2, 1 F2R1][depth - 1]
    alt_hist: Vec<[[[u64; MAX_DEPTH]; 2]; 4]>,
    alt_table: String,
}

impl F1R2Counts {
    pub fn new(sample: &str) -> F1R2Counts {
        F1R2Counts { sample: sample.to_string(), ref_hist: vec![[0; MAX_DEPTH]; 64], alt_hist: vec![[[[0; MAX_DEPTH]; 2]; 4]; 64], alt_table: String::new() }
    }

    /// `F1R2CountsCollector.process` for one locus; `ctx` = reference bases at pos-1, pos, pos+1 (None at the
    /// contig ends).
    pub fn process(&mut self, ctx: Option<[u8; 3]>, c: &LocusCounts) {
        let Some(ctx) = ctx else { return };
        if ctx.contains(&b'N') {
            return;
        }
        let Some(ki) = kmer_index(&ctx) else { return };
        let ref_base = base_index(ctx[1]).unwrap();
        let counts = c.cnt.map(|x| x as i32);
        let depth: i32 = counts.iter().sum();
        // isPileupGood
        let indel_threshold = depth / 100;
        let mut is_indel = c.indel_adj as i32 > indel_threshold;
        is_indel = is_indel || (depth == 0 && c.n_kept > 0);
        if !(depth > 0 && !is_indel) {
            return;
        }
        if !(c.median_mq() >= MIN_MEDIAN_MQ) {
            return;
        }
        let mut copy = counts;
        copy[ref_base] = -1;
        let mut alt = 0;
        for i in 1..4 {
            if copy[i] > copy[alt] {
                alt = i;
            }
        }
        let d = (depth as usize).min(MAX_DEPTH);
        if counts[alt] == 0 {
            self.ref_hist[ki][d - 1] += 1;
            return;
        }
        let (ref_count, alt_count) = (counts[ref_base], counts[alt]);
        let (ref_f1r2, alt_f1r2) = (c.f1r2[ref_base], c.f1r2[alt]);
        if alt_count == 1 {
            let o = if alt_f1r2 == 1 { 0 } else { 1 };
            self.alt_hist[ki][alt][o][d - 1] += 1;
            return;
        }
        use std::fmt::Write as _;
        writeln!(
            self.alt_table,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            std::str::from_utf8(&ctx).unwrap(),
            ref_count,
            alt_count,
            ref_f1r2,
            alt_f1r2,
            ref_count + alt_count,
            BASES[alt] as char
        )
        .unwrap();
    }

    /// Adds the counts of another contig (processed after this one's).
    pub fn merge(&mut self, o: F1R2Counts) {
        for (a, b) in self.ref_hist.iter_mut().zip(&o.ref_hist) {
            for d in 0..MAX_DEPTH {
                a[d] += b[d];
            }
        }
        for (a, b) in self.alt_hist.iter_mut().zip(&o.alt_hist) {
            for x in 0..4 {
                for y in 0..2 {
                    for d in 0..MAX_DEPTH {
                        a[x][y][d] += b[x][y][d];
                    }
                }
            }
        }
        self.alt_table.push_str(&o.alt_table);
    }

    fn metrics_file(&self, labels: &[String], cols: &[&[u64; MAX_DEPTH]]) -> Vec<u8> {
        let mut s = format!("## htsjdk.samtools.metrics.StringHeader\n# {}\n\n\n## HISTOGRAM\tjava.lang.Integer\ndepth", self.sample);
        for l in labels {
            s.push('\t');
            s.push_str(l);
        }
        s.push('\n');
        for d in 0..MAX_DEPTH {
            s.push_str(&(d + 1).to_string());
            for c in cols {
                s.push('\t');
                s.push_str(&c[d].to_string());
            }
            s.push('\n');
        }
        s.push('\n');
        s.into_bytes()
    }

    /// `writeHistograms` + `closeAndArchiveFiles`
    pub fn write_tar_gz(&self, path: &str) {
        let order = kmer_hashmap_order();
        let mut ref_labels = Vec::new();
        let mut ref_cols = Vec::new();
        let mut alt_labels = Vec::new();
        let mut alt_cols = Vec::new();
        for k in &order {
            let ki = kmer_index(k.as_bytes()).unwrap();
            ref_labels.push(k.clone());
            ref_cols.push(&self.ref_hist[ki]);
            let r = ki / 4 % 4;
            for o in [1, 0] {
                for a in 0..4 {
                    if a == r {
                        continue;
                    }
                    alt_labels.push(format!("{k}_{}_{}", BASES[a] as char, if o == 0 { "F1R2" } else { "F2R1" }));
                    alt_cols.push(&self.alt_hist[ki][a][o]);
                }
            }
        }
        let enc = url_encode(&self.sample);
        let table = format!("#<METADATA>SAMPLE={}\ncontext\tref_count\talt_count\tref_f1r2\talt_f1r2\tdepth\talt\n{}", self.sample, self.alt_table);
        write_tar_gz(
            path,
            &[
                (format!("{enc}.ref_histogram"), self.metrics_file(&ref_labels, &ref_cols)),
                (format!("{enc}.alt_table"), table.into_bytes()),
                (format!("{enc}.alt_histogram"), self.metrics_file(&alt_labels, &alt_cols)),
            ],
        );
    }
}

// ------------------------------------------------------------------------------------------------
// LearnReadOrientationModel

/// A metrics file's header line (`# sample`) and histograms (label -> counts by depth 1..=MAX_DEPTH).
fn read_metrics(text: &str) -> (String, Vec<(String, Vec<f64>)>) {
    let lines: Vec<&str> = text.lines().collect();
    let mut sample = String::new();
    let mut i = 0;
    while i < lines.len() {
        let l = lines[i];
        if l.starts_with("## htsjdk.samtools.metrics.StringHeader") && i + 1 < lines.len() {
            sample = lines[i + 1].strip_prefix("# ").unwrap_or(lines[i + 1]).to_string();
        }
        if l.starts_with("## HISTOGRAM") {
            break;
        }
        i += 1;
    }
    let labels: Vec<&str> = lines[i + 1].split('\t').skip(1).collect();
    let mut hists: Vec<(String, Vec<f64>)> = labels.iter().map(|l| (l.to_string(), vec![0.0; MAX_DEPTH])).collect();
    for l in &lines[i + 2..] {
        if l.is_empty() {
            break;
        }
        let f: Vec<&str> = l.split('\t').collect();
        let d: usize = f[0].parse().unwrap();
        if d == 0 || d > MAX_DEPTH {
            continue;
        }
        for (j, v) in f[1..].iter().enumerate() {
            hists[j].1[d - 1] = v.parse().unwrap();
        }
    }
    (sample, hists)
}

#[derive(Clone)]
struct AltRec {
    ctx: String,
    ref_count: i32,
    alt_count: i32,
    ref_f1r2: i32,
    alt_f1r2: i32,
    alt: usize,
}

impl AltRec {
    fn depth(&self) -> i32 {
        self.ref_count + self.alt_count
    }
    fn reverse_complement(&self) -> AltRec {
        AltRec {
            ctx: revcomp(&self.ctx),
            ref_count: self.ref_count,
            alt_count: self.alt_count,
            ref_f1r2: self.ref_count - self.ref_f1r2,
            alt_f1r2: self.alt_count - self.alt_f1r2,
            alt: 3 - self.alt,
        }
    }
}

/// `ArtifactPrior`: pi by state, numExamples, numAltExamples
#[derive(Clone)]
pub struct Prior {
    pub pi: [f64; NUM_STATES],
    pub n: i32,
    pub n_alt: i32,
}

fn learn_context(ctx: &str, ref_hist: &[f64], alt_hists: &[(usize, usize, Vec<f64>)], design: &[AltRec]) -> Prior {
    const CONVERGENCE: f64 = 1e-4;
    const MAX_ITER: usize = 20;
    let ref_allele = base_index(ctx.as_bytes()[1]).unwrap();
    let num_alt = design.len() as i32 + alt_hists.iter().map(|h| h.2.iter().sum::<f64>() as i32).sum::<i32>();
    let num_ref = ref_hist.iter().sum::<f64>() as i32;
    let pseudo = flat_prior(ref_allele);
    let mut prior = pseudo;
    let mut iters = 0;
    let mut ref_resp = vec![[0.0; NUM_STATES]; MAX_DEPTH];
    let mut alt_resp = vec![[0.0; NUM_STATES]; design.len()];
    let mut d1_resp = vec![[[[0.0; NUM_STATES]; 2]; 4]; MAX_DEPTH];
    loop {
        let old = prior;
        // E step
        for i in 0..MAX_DEPTH {
            ref_resp[i] = compute_responsibilities(ref_allele, ref_allele, 0, 0, (i + 1) as i32, &prior, false);
        }
        for (n, r) in design.iter().enumerate() {
            alt_resp[n] = compute_responsibilities(ref_allele, r.alt, r.alt_count, r.alt_f1r2, r.depth(), &prior, false);
        }
        for i in 0..MAX_DEPTH {
            for a in 0..4 {
                if a == ref_allele {
                    continue;
                }
                for o in 0..2 {
                    d1_resp[i][a][o] = compute_responsibilities(ref_allele, a, 1, if o == 0 { 1 } else { 0 }, (i + 1) as i32, &prior, false);
                }
            }
        }
        // M step
        let mut eff_design = [0.0; NUM_STATES];
        if !design.is_empty() {
            eff_design = alt_resp[0];
            for row in &alt_resp[1..] {
                for s in 0..NUM_STATES {
                    eff_design[s] += row[s];
                }
            }
        }
        let mut eff_hist = [0.0; NUM_STATES];
        for (a, o, h) in alt_hists {
            let mut e = [0.0; NUM_STATES];
            for i in 0..MAX_DEPTH {
                for s in 0..NUM_STATES {
                    let v = h[i] * d1_resp[i][*a][*o][s];
                    if i == 0 {
                        e[s] = v;
                    } else {
                        e[s] += v;
                    }
                }
            }
            for s in 0..NUM_STATES {
                eff_hist[s] += e[s];
            }
        }
        let mut eff_ref = [0.0; NUM_STATES];
        for i in 0..MAX_DEPTH {
            for s in 0..NUM_STATES {
                let v = ref_hist[i] * ref_resp[i][s];
                if i == 0 {
                    eff_ref[s] = v;
                } else {
                    eff_ref[s] += v;
                }
            }
        }
        let mut c = [0.0; NUM_STATES];
        for s in 0..NUM_STATES {
            c[s] = ((eff_design[s] + eff_hist[s]) + eff_ref[s]) + pseudo[s];
        }
        let mut sum = 0.0;
        for &x in &c {
            sum += x;
        }
        for s in 0..NUM_STATES {
            prior[s] = c[s] / sum;
        }
        let mut d2 = 0.0;
        for s in 0..NUM_STATES {
            let dp = old[s] - prior[s];
            d2 += dp * dp;
        }
        iters += 1;
        if !(d2.sqrt() > CONVERGENCE && iters < MAX_ITER) {
            break;
        }
    }
    Prior { pi: prior, n: num_alt + num_ref, n_alt: num_alt }
}

fn write_priors(sample: &str, priors: &std::collections::HashMap<String, Prior>) -> Vec<u8> {
    let mut s = format!(
        "#<METADATA>SAMPLE={sample}\ncontext\trev_comp\tf1r2_a\tf1r2_c\tf1r2_g\tf1r2_t\tf2r1_a\tf2r1_c\tf2r1_g\tf2r1_t\thom_ref\tgermline_het\tsomatic_het\thom_var\tnum_examples\tnum_alt_examples\n"
    );
    for k in kmer_hashmap_order() {
        let p = &priors[&k];
        s.push_str(&k);
        s.push('\t');
        s.push_str(&revcomp(&k));
        for x in p.pi {
            s.push('\t');
            s.push_str(&java_double_to_string(x));
        }
        s.push_str(&format!("\t{}\t{}\n", p.n, p.n_alt));
    }
    s.into_bytes()
}

/// `ArtifactPriorCollection.set`: the canonical prior and its reverse complement
fn set_prior(map: &mut std::collections::HashMap<String, Prior>, ctx: &str, p: Prior) {
    let mut rc = p.clone();
    for s in 0..NUM_STATES {
        rc.pi[s] = p.pi[rev_comp_state(s)];
    }
    map.insert(ctx.to_string(), p);
    map.insert(revcomp(ctx), rc);
}

fn flat_collection() -> std::collections::HashMap<String, Prior> {
    all_kmers().into_iter().map(|k| {
        let r = base_index(k.as_bytes()[1]).unwrap();
        (k, Prior { pi: flat_prior(r), n: 0, n_alt: 0 })
    }).collect()
}

const LROM_USAGE: &str = "usage: variant-rs learnreadorientationmodel -I f1r2.tar.gz [-I ...] -O artifact-priors.tar.gz";

pub fn run_lrom(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let (mut inputs, mut output) = (Vec::<String>::new(), None::<String>);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{LROM_USAGE}"));
        match argv[i].as_str() {
            "-I" | "--input" => inputs.push(v()),
            "-O" | "--output" => output = Some(v()),
            a => panic!("unknown argument {a}\n{LROM_USAGE}"),
        }
        i += 2;
    }
    let output = output.expect(LROM_USAGE);
    assert!(output.ends_with(".tar.gz"), "Output file must end in .tar.gz");
    // per sample: ref histograms, alt histograms (summed over files), alt records (concatenated)
    type Hists = Vec<(String, Vec<f64>)>;
    let mut ref_h: Vec<(String, Hists)> = Vec::new();
    let mut alt_h: Vec<(String, Hists)> = Vec::new();
    let mut recs: Vec<(String, Vec<AltRec>)> = Vec::new();
    let add = |dst: &mut Vec<(String, Hists)>, sample: String, h: Hists| match dst.iter_mut().find(|(s, _)| *s == sample) {
        None => dst.push((sample, h)),
        Some((_, acc)) => {
            for (l, v) in h {
                let t = acc.iter_mut().find(|(x, _)| *x == l).unwrap_or_else(|| panic!("Missing histogram header for: {l}"));
                for d in 0..MAX_DEPTH {
                    t.1[d] += v[d];
                }
            }
        }
    };
    for path in &inputs {
        for (name, data) in read_tar_gz(path) {
            let text = String::from_utf8_lossy(&data).to_string();
            if name.ends_with(".ref_histogram") {
                let (s, h) = read_metrics(&text);
                add(&mut ref_h, s, h);
            } else if name.ends_with(".alt_histogram") {
                let (s, h) = read_metrics(&text);
                add(&mut alt_h, s, h);
            } else if name.ends_with(".alt_table") {
                let mut sample = String::new();
                let mut v = Vec::new();
                let mut header_seen = false;
                for l in text.lines() {
                    if let Some(m) = l.strip_prefix("#<METADATA>SAMPLE=") {
                        sample = m.to_string();
                        continue;
                    }
                    if l.starts_with('#') || l.is_empty() {
                        continue;
                    }
                    if !header_seen {
                        header_seen = true;
                        continue;
                    }
                    let f: Vec<&str> = l.split('\t').collect();
                    v.push(AltRec {
                        ctx: f[0].to_string(),
                        ref_count: f[1].parse().unwrap(),
                        alt_count: f[2].parse().unwrap(),
                        ref_f1r2: f[3].parse().unwrap(),
                        alt_f1r2: f[4].parse().unwrap(),
                        alt: base_index(f[6].as_bytes()[0]).unwrap(),
                    });
                }
                match recs.iter_mut().find(|(s, _)| *s == sample) {
                    Some((_, acc)) => acc.extend(v),
                    None => recs.push((sample, v)),
                }
            }
        }
    }
    let mut files = Vec::new();
    for (sample, records) in &recs {
        let refs = &ref_h.iter().find(|(s, _)| s == sample).expect("no ref histograms for sample").1;
        let alts = alt_h.iter().find(|(s, _)| s == sample).map(|x| &x.1);
        let mut coll = flat_collection();
        let mut n_learned = 0;
        for ctx in canonical_kmers() {
            let rc = revcomp(&ctx);
            let get_ref = |k: &str| refs.iter().find(|(l, _)| l == k).map_or(vec![0.0; MAX_DEPTH], |x| x.1.clone());
            let (r1, r2) = (get_ref(&ctx), get_ref(&rc));
            let ref_hist: Vec<f64> = (0..MAX_DEPTH).map(|d| r1[d] + r2[d]).collect();
            // combineAltDepthOneHistogramWithRC
            let mut alt_hists: Vec<(usize, usize, Vec<f64>)> = Vec::new();
            if let Some(alts) = alts {
                let ra = base_index(ctx.as_bytes()[1]).unwrap();
                for a in 0..4 {
                    if a == ra {
                        continue;
                    }
                    for o in 0..2 {
                        let on = ["F1R2", "F2R1"];
                        let l1 = format!("{ctx}_{}_{}", BASES[a] as char, on[o]);
                        let l2 = format!("{rc}_{}_{}", BASES[3 - a] as char, on[1 - o]);
                        let g = |l: &str| alts.iter().find(|(x, _)| x == l).map_or(vec![0.0; MAX_DEPTH], |x| x.1.clone());
                        let (h1, h2) = (g(&l1), g(&l2));
                        alt_hists.push((a, o, (0..MAX_DEPTH).map(|d| h1[d] + h2[d]).collect()));
                    }
                }
            }
            let mut design: Vec<AltRec> = records.iter().filter(|r| r.ctx == ctx).cloned().collect();
            design.extend(records.iter().filter(|r| r.ctx == rc).map(|r| r.reverse_complement()));
            if ref_hist.iter().sum::<f64>() == 0.0 || design.is_empty() {
                continue;
            }
            let p = learn_context(&ctx, &ref_hist, &alt_hists, &design);
            set_prior(&mut coll, &ctx, p);
            n_learned += 1;
        }
        eprintln!("learnreadorientationmodel: {sample}: {n_learned} contexts learned");
        files.push((format!("{}.orientation_priors", url_encode(sample)), write_priors(sample, &coll)));
    }
    write_tar_gz(&output, &files);
    eprintln!("learnreadorientationmodel: {:.2}s", t0.elapsed().as_secs_f64());
}

/// `ArtifactPriorCollection.readArtifactPriors` of every table in the .tar.gz files: (sample, priors by k-mer)
pub fn read_priors(paths: &[String]) -> Vec<(String, std::collections::HashMap<String, Prior>)> {
    let mut out = Vec::new();
    for path in paths {
        for (_, data) in read_tar_gz(path) {
            let text = String::from_utf8_lossy(&data).to_string();
            let mut sample = String::new();
            let mut header: Option<Vec<String>> = None;
            let mut rows: Vec<(String, Prior)> = Vec::new();
            for l in text.lines() {
                if let Some(m) = l.strip_prefix("#<METADATA>SAMPLE=") {
                    sample = m.to_string();
                    continue;
                }
                if l.starts_with('#') || l.is_empty() {
                    continue;
                }
                let f: Vec<&str> = l.split('\t').collect();
                match &header {
                    None => header = Some(f.iter().map(|s| s.to_string()).collect()),
                    Some(h) => {
                        let col = |n: &str| f[h.iter().position(|x| x == n).unwrap()];
                        let names = ["f1r2_a", "f1r2_c", "f1r2_g", "f1r2_t", "f2r1_a", "f2r1_c", "f2r1_g", "f2r1_t", "hom_ref", "germline_het", "somatic_het", "hom_var"];
                        let mut pi = [0.0; NUM_STATES];
                        for (s, n) in names.iter().enumerate() {
                            pi[s] = col(n).parse().unwrap();
                        }
                        rows.push((col("context").to_string(), Prior { pi, n: col("num_examples").parse().unwrap(), n_alt: col("num_alt_examples").parse().unwrap() }));
                    }
                }
            }
            let mut coll = flat_collection();
            for ctx in canonical_kmers() {
                let p = rows.iter().find(|(c, _)| *c == ctx).unwrap_or_else(|| panic!("ArtifactPrior object isn't present for reference context {ctx}")).1.clone();
                set_prior(&mut coll, &ctx, p);
            }
            out.push((sample, coll));
        }
    }
    out
}

#[allow(dead_code)]
fn _kmer_of_unused() -> String {
    kmer_of(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn kmers() {
        let k = all_kmers();
        assert_eq!(k[0], "AAC");
        assert_eq!(k[63], "AAA");
        assert_eq!(canonical_kmers().len(), 32);
        assert_eq!(rev_comp_state(0), 7);
        assert_eq!(rev_comp_state(5), 2);
        assert_eq!(url_encode("tumor sample"), "tumor+sample");
    }
}
