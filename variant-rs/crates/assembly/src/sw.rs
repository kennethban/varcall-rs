//! Port of GATK's SmithWatermanJavaAligner (which the GKL AVX aligner reproduces),
//! including the exact-substring shortcut of SWNativeAlignerWrapper / SmithWatermanJavaAligner.
//!
//! The matrices are filled by anti-diagonals, as GKL's `smithWatermanBackTrack` does: the cells of one
//! anti-diagonal depend only on the two previous ones, so the inner loop has no loop-carried dependency and
//! the compiler vectorises it. Every cell evaluates the Java recurrence with the same operations and tie
//! rules, so the scores are those of the Java row-by-row fill. Like GKL, the traceback matrix holds one byte
//! per cell: the winning move and whether each gap state extended the previous cell's gap; the Java gap sizes
//! are the lengths of those extension chains. The end-cell selection and the traceback are the Java ones.

use crate::cigar::Cigar;

#[derive(Clone, Copy, Debug)]
pub struct SwParameters {
    pub match_value: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

/// SmithWatermanAlignmentConstants.STANDARD_NGS (dangling end default)
pub const STANDARD_NGS: SwParameters =
    SwParameters { match_value: 25, mismatch_penalty: -50, gap_open_penalty: -110, gap_extend_penalty: -6 };
/// SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS (haplotype to reference default)
pub const NEW_SW_PARAMETERS: SwParameters =
    SwParameters { match_value: 200, mismatch_penalty: -150, gap_open_penalty: -260, gap_extend_penalty: -11 };

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Overhang {
    SoftClip,
    Indel,
    LeadingIndel,
    Ignore,
}

pub struct SwAlignment {
    pub cigar: Cigar,
    pub offset: i32,
}

/// The column a fill resumes after: scores H(0..=n, c) and horizontal gaps E(0..=n, c).
struct Resume<'a> {
    c: usize,
    h: &'a [i32],
    e: &'a [i32],
}

/// Reusable work buffers.
#[derive(Default)]
pub struct SwAligner {
    /// traceback byte of cell (i, j) at `(i + j) * (n + 1) + i` (anti-diagonal major)
    bt: Vec<u8>,
    /// scores of the last three anti-diagonals, indexed by row
    h: [Vec<i32>; 3],
    /// best horizontal gap ending in the last cell computed in each row (Java `bestGapH`)
    e: Vec<i32>,
    /// best vertical gap ending in the last cell computed in each column (Java `bestGapV`), column j at m - j
    f: Vec<i32>,
    /// the alternate reversed: column j's base at m - j
    alt_rev: Vec<u8>,
    /// scores of the last column (j = m) and the last row (i = n), all `calculateCigar` reads
    last_col: Vec<i32>,
    last_row: Vec<i32>,
}

/// Utils.lastIndexOf
pub fn last_index_of(reference: &[u8], query: &[u8]) -> i32 {
    let ql = query.len();
    if ql > reference.len() {
        return -1;
    }
    let mut r = reference.len() - ql;
    loop {
        if &reference[r..r + ql] == query {
            return r as i32;
        }
        if r == 0 {
            return -1;
        }
        r -= 1;
    }
}

/// `VARIANT_RS_SW_DUMP=<file>`: records every aligner call (for the `swbench` replay benchmark). Lines:
/// `S <params> <strategy> <ref> <alt>` (SwAligner), and per FixedRefAligner reference a block
/// `R <params> <ref>`, `B` (start of one graph's haplotype batch), `A <alt>`.
pub mod dump {
    use std::io::Write;
    use std::sync::{Mutex, OnceLock};
    fn w() -> &'static Option<Mutex<std::io::BufWriter<std::fs::File>>> {
        static W: OnceLock<Option<Mutex<std::io::BufWriter<std::fs::File>>>> = OnceLock::new();
        W.get_or_init(|| {
            std::env::var_os("VARIANT_RS_SW_DUMP").map(|p| Mutex::new(std::io::BufWriter::new(std::fs::File::create(p).unwrap())))
        })
    }
    #[inline]
    pub fn on() -> bool {
        w().is_some()
    }
    pub fn write(s: &[u8]) {
        if let Some(m) = w() {
            m.lock().unwrap().write_all(s).unwrap();
        }
    }
    pub fn flush() {
        if let Some(m) = w() {
            m.lock().unwrap().flush().unwrap();
        }
    }
    pub fn params(p: &super::SwParameters) -> String {
        format!("{},{},{},{}", p.match_value, p.mismatch_penalty, p.gap_open_penalty, p.gap_extend_penalty)
    }
}

const MATRIX_MIN_CUTOFF: i32 = -100_000_000;
const LOW_INIT_VALUE: i32 = i32::MIN / 2;

// traceback byte: winning move in the low two bits, gap-extension flags above
const INS: u8 = 1; // from the horizontal gap (Java step right, a negative btrack value)
const DEL: u8 = 2; // from the vertical gap (Java step down, a positive btrack value)
const E_EXT: u8 = 4; // the horizontal gap ending here extends the one ending in (i, j - 1)
const F_EXT: u8 = 8; // the vertical gap ending here extends the one ending in (i - 1, j)

/// One anti-diagonal: cells `k` of every slice are the same cell (i, j), with i increasing and j decreasing.
/// `diag`, `top` and `left` are the scores of (i-1, j-1), (i-1, j) and (i, j-1); `e` and `f` hold the gaps
/// ending in (i, j-1) and (i-1, j) and receive those ending in (i, j).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn antidiagonal(
    diag: &[i32],
    top: &[i32],
    left: &[i32],
    out: &mut [i32],
    e: &mut [i32],
    f: &mut [i32],
    ref_bases: &[u8],
    alt_bases: &[u8],
    bt: &mut [u8],
    p: &SwParameters,
) {
    let len = out.len();
    let (diag, top, left, e, f, ref_bases, alt_bases, bt) =
        (&diag[..len], &top[..len], &left[..len], &mut e[..len], &mut f[..len], &ref_bases[..len], &alt_bases[..len], &mut bt[..len]);
    let (w_match, w_mismatch, w_open, w_extend) = (p.match_value, p.mismatch_penalty, p.gap_open_penalty, p.gap_extend_penalty);
    for k in 0..len {
        let step_diag = diag[k] + if ref_bases[k] == alt_bases[k] { w_match } else { w_mismatch };
        // vertical gap (step down)
        let open_v = top[k] + w_open;
        let ext_v = f[k] + w_extend;
        let (step_down, fx) = if open_v > ext_v { (open_v, 0) } else { (ext_v, F_EXT) };
        f[k] = step_down;
        // horizontal gap (step right)
        let open_h = left[k] + w_open;
        let ext_h = e[k] + w_extend;
        let (step_right, ex) = if open_h > ext_h { (open_h, 0) } else { (ext_h, E_EXT) };
        e[k] = step_right;
        let (v, mv) = if step_diag >= step_down && step_diag >= step_right {
            (step_diag, 0)
        } else if step_right >= step_down {
            (step_right, INS)
        } else {
            (step_down, DEL)
        };
        out[k] = MATRIX_MIN_CUTOFF.max(v);
        bt[k] = mv | ex | fx;
    }
}

impl SwAligner {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn align(&mut self, reference: &[u8], alternate: &[u8], p: &SwParameters, strategy: Overhang) -> SwAlignment {
        assert!(!reference.is_empty() && !alternate.is_empty(), "Non-null, non-empty sequences are required");
        let _p = crate::prof::scope("SwAligner.align (all callers)");
        if dump::on() {
            let mut l = format!("S\t{}\t{:?}\t", dump::params(p), strategy).into_bytes();
            l.extend_from_slice(reference);
            l.push(b'\t');
            l.extend_from_slice(alternate);
            l.push(b'\n');
            dump::write(&l);
        }
        if strategy == Overhang::SoftClip || strategy == Overhang::Ignore {
            let idx = last_index_of(reference, alternate);
            if idx != -1 {
                return SwAlignment { cigar: vec![(b'M', alternate.len() as u32)], offset: idx };
            }
        }
        crate::prof::add("SwAligner calls (DP)", 1);
        crate::prof::add("SwAligner cells", (reference.len() * alternate.len()) as u64);
        self.align_dp(reference, alternate, p, strategy)
    }

    /// Smith-Waterman without the exact-substring shortcut.
    fn align_dp(&mut self, reference: &[u8], alternate: &[u8], p: &SwParameters, strategy: Overhang) -> SwAlignment {
        self.fill(reference, alternate, p, strategy, None, &[], &mut Vec::new());
        self.calculate_cigar(reference.len(), alternate.len(), strategy)
    }

    /// Fills the columns after `resume.c` (all columns without `resume`). The traceback cells and last-row scores
    /// of the columns up to `c` are kept from earlier fills, which must have had the same reference and the same
    /// first `c` alternate bases. For every column `t` in `save` (all > c), writes its scores H(0..=n, t) and
    /// then its horizontal gaps E(0..=n, t) to `saved`, `2 (n + 1)` values per column, for resuming later.
    #[allow(clippy::too_many_arguments)]
    fn fill(
        &mut self,
        reference: &[u8],
        alternate: &[u8],
        p: &SwParameters,
        strategy: Overhang,
        resume: Option<Resume>,
        save: &[usize],
        saved: &mut Vec<i32>,
    ) {
        let n = reference.len();
        let m = alternate.len();
        let stride = n + 1;
        let c = resume.as_ref().map_or(0, |r| r.c);
        assert!(c <= m && save.iter().all(|&t| t > c && t <= m));
        // the sub-matrix of columns c+1..=m, with local column j - c and local anti-diagonal d = i + j - c
        let mw = m - c;
        if self.bt.len() < (n + m + 1) * stride {
            self.bt.resize((n + m + 1) * stride, 0);
        }
        for h in self.h.iter_mut() {
            h.clear();
            h.resize(n + 1, 0);
        }
        self.e.clear();
        match &resume {
            Some(r) => self.e.extend_from_slice(&r.e[..=n]),
            None => self.e.resize(n + 1, LOW_INIT_VALUE),
        }
        self.f.clear();
        self.f.resize(mw, LOW_INIT_VALUE);
        self.alt_rev.clear();
        self.alt_rev.extend(alternate[c..].iter().rev());
        self.last_col.clear();
        self.last_col.resize(n + 1, 0);
        self.last_row.resize(m + 1, 0);
        saved.clear();
        saved.resize(save.len() * 2 * stride, 0);
        // row 0 and column 0: 0, or gap penalties for the INDEL strategies; the left column when resuming
        let indel = strategy == Overhang::Indel || strategy == Overhang::LeadingIndel;
        let border = |k: usize| if indel && k > 0 { p.gap_open_penalty + (k as i32 - 1) * p.gap_extend_penalty } else { 0 };
        let left = |i: usize| match &resume {
            Some(r) => r.h[i],
            None => border(i),
        };
        for (q, &t) in save.iter().enumerate() {
            saved[q * 2 * stride] = border(t);
        }
        if mw == 0 {
            for i in 1..=n {
                self.last_col[i] = left(i);
            }
            return;
        }
        let [h2, h1, h0] = &mut self.h;
        let (mut h2, mut h1, mut h0) = (std::mem::take(h2), std::mem::take(h1), std::mem::take(h0));
        // local anti-diagonals 0 and 1 are border cells: (0,c); (0,c+1) and (1,c)
        h2[0] = border(c);
        h1[0] = border(c + 1);
        h1[1] = left(1);
        for d in 2..=n + mw {
            let lo = if d > mw { d - mw } else { 1 };
            let hi = n.min(d - 1);
            let len = hi + 1 - lo;
            let fo = mw + lo - d; // index of local column d - lo in `f` and `alt_rev`
            let bo = (d + c) * stride + lo;
            antidiagonal(
                &h2[lo - 1..hi],
                &h1[lo - 1..hi],
                &h1[lo..=hi],
                &mut h0[lo..=hi],
                &mut self.e[lo..=hi],
                &mut self.f[fo..fo + len],
                &reference[lo - 1..hi],
                &self.alt_rev[fo..fo + len],
                &mut self.bt[bo..bo + len],
                p,
            );
            h0[0] = border(c + d);
            if d <= n {
                h0[d] = left(d);
            }
            if d > mw {
                self.last_col[d - mw] = h0[d - mw];
            }
            if d > n {
                self.last_row[c + d - n] = h0[n];
            }
            for (q, &t) in save.iter().enumerate() {
                let tl = t - c;
                if d > tl && d - tl <= n {
                    let i = d - tl;
                    saved[q * 2 * stride + i] = h0[i];
                    saved[q * 2 * stride + stride + i] = self.e[i];
                }
            }
            (h2, h1, h0) = (h1, h0, h2);
        }
        self.h = [h2, h1, h0];
    }

    fn calculate_cigar(&self, ref_len: usize, alt_len: usize, strategy: Overhang) -> SwAlignment {
        let ref_length = ref_len as i64;
        let alt_length = alt_len as i64;
        let mut p1: i64 = 0;
        let mut p2: i64;
        let mut maxscore = i32::MIN;
        let mut segment_length: i64 = 0;
        if strategy == Overhang::Indel {
            p1 = ref_length;
            p2 = alt_length;
        } else {
            p2 = alt_length;
            for i in 1..=ref_len {
                let cur = self.last_col[i];
                if cur >= maxscore {
                    p1 = i as i64;
                    maxscore = cur;
                }
            }
            if strategy != Overhang::LeadingIndel {
                for j in 1..=alt_len {
                    let cur = self.last_row[j];
                    if cur > maxscore || (cur == maxscore && (ref_length - j as i64).abs() < (p1 - p2).abs()) {
                        p1 = ref_length;
                        p2 = j as i64;
                        maxscore = cur;
                        segment_length = alt_length - j as i64;
                    }
                }
            }
        }
        self.traceback(ref_len + 1, p1, p2, segment_length, strategy)
    }

    /// The Java btrack value of cell (i, j): 0 for a match, the vertical gap size for a deletion, minus the
    /// horizontal gap size for an insertion.
    fn btrack(&self, stride: usize, i: usize, j: usize) -> i32 {
        let at = |i: usize, j: usize| self.bt[(i + j) * stride + i];
        let b = at(i, j);
        if b & 3 == INS {
            let (mut k, mut c) = (1, j);
            while at(i, c) & E_EXT != 0 {
                k += 1;
                c -= 1;
            }
            -k
        } else if b & 3 == DEL {
            let (mut k, mut r) = (1, i);
            while at(r, j) & F_EXT != 0 {
                k += 1;
                r -= 1;
            }
            k
        } else {
            0
        }
    }

    /// The traceback of `SmithWatermanJavaAligner.calculateCigar` from the end cell `(p1, p2)`.
    fn traceback(&self, stride: usize, mut p1: i64, mut p2: i64, mut segment_length: i64, strategy: Overhang) -> SwAlignment {
        let mut lce: Vec<(u8, u32)> = Vec::with_capacity(8);
        if segment_length > 0 && strategy == Overhang::SoftClip {
            lce.push((b'S', segment_length as u32));
            segment_length = 0;
        }
        // states: M=match, I=insertion, D=deletion
        let mut state = b'M';
        loop {
            let btr = self.btrack(stride, p1 as usize, p2 as usize);
            let new_state;
            let mut step_length: i64 = 1;
            if btr > 0 {
                new_state = b'D';
                step_length = btr as i64;
            } else if btr < 0 {
                new_state = b'I';
                step_length = -(btr as i64);
            } else {
                new_state = b'M';
            }
            match new_state {
                b'M' => {
                    p1 -= 1;
                    p2 -= 1;
                }
                b'I' => p2 -= step_length,
                _ => p1 -= step_length,
            }
            if new_state == state {
                segment_length += step_length;
            } else {
                if segment_length > 0 {
                    lce.push((state, segment_length as u32));
                }
                segment_length = step_length;
                state = new_state;
            }
            if !(p1 > 0 && p2 > 0) {
                break;
            }
        }
        let offset;
        match strategy {
            Overhang::SoftClip => {
                lce.push((state, segment_length as u32));
                if p2 > 0 {
                    lce.push((b'S', p2 as u32));
                }
                offset = p1 as i32;
            }
            Overhang::Ignore => {
                lce.push((state, (segment_length + p2) as u32));
                offset = (p1 - p2) as i32;
            }
            _ => {
                lce.push((state, segment_length as u32));
                if p1 > 0 {
                    lce.push((b'D', p1 as u32));
                } else if p2 > 0 {
                    lce.push((b'I', p2 as u32));
                }
                offset = 0;
            }
        }
        lce.reverse();
        // htsjdk Cigar keeps elements as given (zero-length elements possible but not produced here)
        SwAlignment { cigar: lce, offset }
    }
}

/// Haplotype-to-reference alignments (SOFTCLIP strategy) against one fixed reference.
///
/// Column j of the matrices depends only on the first j alternate bases, so a fill resumes after the prefix an
/// alternate shares with the previously aligned one: the traceback cells and last-row scores of those columns
/// are still in place (cells sit at fixed positions for a fixed reference), and the left boundary comes from a
/// saved column (scores and horizontal gaps of every row). Which columns to save is known from the rest of
/// the batch: a column is needed while later alternates still share its prefix, so the saved columns form a
/// stack (a depth-first walk of the alternates' trie when the batch is sorted). Results are those of aligning
/// each alternate from scratch.
#[derive(Default)]
pub struct FixedRefAligner {
    reference: Vec<u8>,
    params: Option<SwParameters>,
    sw: SwAligner,
    /// the last alternate filled: traceback columns 1..=alt.len() and the last-row scores belong to it
    alt: Vec<u8>,
    /// saved columns (increasing) and their data, `2 (n + 1)` values each (see `SwAligner::fill`)
    ck_cols: Vec<usize>,
    ck: Vec<i32>,
    saved: Vec<i32>,
    /// VARIANT_RS_SW_DUMP block for the current reference
    log: Vec<u8>,
}

impl Drop for FixedRefAligner {
    fn drop(&mut self) {
        self.flush_log();
    }
}

fn common_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}

impl FixedRefAligner {
    pub fn new() -> Self {
        Self::default()
    }

    fn flush_log(&mut self) {
        if !self.log.is_empty() {
            dump::write(&self.log);
            self.log.clear();
        }
    }

    /// Marks the start of one graph's batch of haplotypes (VARIANT_RS_SW_DUMP only).
    pub fn begin_batch(&mut self) {
        if dump::on() {
            self.log.extend_from_slice(b"B\n");
        }
    }

    pub fn set_reference(&mut self, reference: &[u8], params: &SwParameters) {
        if dump::on() {
            self.flush_log();
            self.log.extend_from_slice(format!("R\t{}\t", dump::params(params)).as_bytes());
            self.log.extend_from_slice(reference);
            self.log.push(b'\n');
        }
        self.reference.clear();
        self.reference.extend_from_slice(reference);
        self.params = Some(*params);
        self.alt.clear();
        self.ck_cols.clear();
        self.ck.clear();
    }

    pub fn align(&mut self, alternate: &[u8]) -> SwAlignment {
        self.align_batch(&[alternate]).pop().unwrap()
    }

    /// Aligns each alternate (as `align` would, one after the other). Sorted batches share the most columns.
    pub fn align_batch(&mut self, alts: &[&[u8]]) -> Vec<SwAlignment> {
        assert!(self.params.is_some(), "reference not set");
        assert!(!self.reference.is_empty());
        let mut out: Vec<Option<SwAlignment>> = Vec::with_capacity(alts.len());
        let mut dp: Vec<usize> = Vec::new();
        for (k, a) in alts.iter().enumerate() {
            assert!(!a.is_empty());
            if dump::on() {
                self.log.extend_from_slice(b"A\t");
                self.log.extend_from_slice(a);
                self.log.push(b'\n');
            }
            let idx = last_index_of(&self.reference, a);
            if idx != -1 {
                out.push(Some(SwAlignment { cigar: vec![(b'M', a.len() as u32)], offset: idx }));
            } else {
                out.push(None);
                dp.push(k);
            }
        }
        // prefix each DP alternate shares with the next one
        let next_lcp: Vec<usize> = dp.windows(2).map(|w| common_prefix(alts[w[0]], alts[w[1]])).collect();
        for (q, &k) in dp.iter().enumerate() {
            out[k] = Some(self.align_dp(alts[k], &next_lcp[q.min(next_lcp.len())..]));
        }
        out.into_iter().map(|a| a.unwrap()).collect()
    }

    /// Smith-Waterman of `alternate`. `future` lists the prefixes shared by consecutive DP alternates from this
    /// one to the end of the batch (`future[0]`: this one and the next).
    fn align_dp(&mut self, alternate: &[u8], future: &[usize]) -> SwAlignment {
        let p = self.params.unwrap();
        let n = self.reference.len();
        let m = alternate.len();
        let w = 2 * (n + 1);
        // resume after the deepest saved column within the prefix shared with the last alternate
        let lcp = common_prefix(&self.alt, alternate);
        while self.ck_cols.last().is_some_and(|&c| c > lcp) {
            self.ck_cols.pop();
            self.ck.truncate(self.ck_cols.len() * w);
        }
        let c = self.ck_cols.last().copied().unwrap_or(0);
        // columns to save: the prefixes later alternates share with this one (running minima of `future`)
        // beyond c, increasing
        let mut save: Vec<usize> = Vec::new();
        let mut r = usize::MAX;
        for &t in future {
            if t < r {
                r = t;
                if r <= c {
                    break;
                }
                save.push(r);
            }
        }
        save.reverse();
        crate::prof::add("hapSW calls (DP)", 1);
        crate::prof::add("hapSW cells", ((m - c) * n) as u64);
        crate::prof::add("hapSW cells without prefix reuse", (m * n) as u64);
        let resume = if c > 0 {
            let at = self.ck.len() - w;
            Some(Resume { c, h: &self.ck[at..at + n + 1], e: &self.ck[at + n + 1..at + w] })
        } else {
            None
        };
        self.sw.fill(&self.reference, alternate, &p, Overhang::SoftClip, resume, &save, &mut self.saved);
        for (q, &t) in save.iter().enumerate() {
            self.ck_cols.push(t);
            self.ck.extend_from_slice(&self.saved[q * w..(q + 1) * w]);
        }
        self.alt.clear();
        self.alt.extend_from_slice(alternate);
        self.sw.calculate_cigar(n, m, Overhang::SoftClip)
    }
}
