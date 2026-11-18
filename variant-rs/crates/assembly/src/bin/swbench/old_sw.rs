//! Port of GATK's SmithWatermanJavaAligner (which the GKL AVX aligner reproduces),
//! including the exact-substring shortcut of SWNativeAlignerWrapper / SmithWatermanJavaAligner.

type Cigar = Vec<(u8, u32)>;

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

/// Reusable work buffers.
#[derive(Default)]
pub struct SwAligner {
    sw: Vec<i32>,
    btrack: Vec<i32>,
    best_gap_v: Vec<i32>,
    gap_size_v: Vec<i32>,
    best_gap_h: Vec<i32>,
    gap_size_h: Vec<i32>,
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

const MATRIX_MIN_CUTOFF: i32 = -100_000_000;
const LOW_INIT_VALUE: i32 = i32::MIN / 2;

impl SwAligner {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn align(&mut self, reference: &[u8], alternate: &[u8], p: &SwParameters, strategy: Overhang) -> SwAlignment {
        assert!(!reference.is_empty() && !alternate.is_empty(), "Non-null, non-empty sequences are required");
        if strategy == Overhang::SoftClip || strategy == Overhang::Ignore {
            let idx = last_index_of(reference, alternate);
            if idx != -1 {
                return SwAlignment { cigar: vec![(b'M', alternate.len() as u32)], offset: idx };
            }
        }
        self.calculate_matrix(reference, alternate, p, strategy);
        self.calculate_cigar(reference.len(), alternate.len(), strategy)
    }

    fn calculate_matrix(&mut self, reference: &[u8], alternate: &[u8], p: &SwParameters, strategy: Overhang) {
        let nrow = reference.len() + 1;
        let ncol = alternate.len() + 1;
        let total = nrow * ncol;
        self.sw.clear();
        self.sw.resize(total, 0);
        self.btrack.clear();
        self.btrack.resize(total, 0);
        self.best_gap_v.clear();
        self.best_gap_v.resize(ncol + 1, LOW_INIT_VALUE);
        self.gap_size_v.clear();
        self.gap_size_v.resize(ncol + 1, 0);
        self.best_gap_h.clear();
        self.best_gap_h.resize(nrow + 1, LOW_INIT_VALUE);
        self.gap_size_h.clear();
        self.gap_size_h.resize(nrow + 1, 0);

        let sw = &mut self.sw[..];
        if strategy == Overhang::Indel || strategy == Overhang::LeadingIndel {
            sw[1] = p.gap_open_penalty;
            let mut cur = p.gap_open_penalty;
            for i in 2..ncol {
                cur += p.gap_extend_penalty;
                sw[i] = cur;
            }
            sw[ncol] = p.gap_open_penalty;
            cur = p.gap_open_penalty;
            for i in 2..nrow {
                cur += p.gap_extend_penalty;
                sw[i * ncol] = cur;
            }
        }
        let w_open = p.gap_open_penalty;
        let w_extend = p.gap_extend_penalty;
        let w_match = p.match_value;
        let w_mismatch = p.mismatch_penalty;
        let best_gap_v = &mut self.best_gap_v[..];
        let gap_size_v = &mut self.gap_size_v[..];
        let btrack = &mut self.btrack[..];
        for i in 1..nrow {
            let a_base = reference[i - 1];
            let (prev_rows, cur_rows) = sw.split_at_mut(i * ncol);
            let last_row = &prev_rows[(i - 1) * ncol..];
            let cur_row = &mut cur_rows[..ncol];
            let bt_row = &mut btrack[i * ncol..(i + 1) * ncol];
            let mut best_gap_h = self.best_gap_h[i];
            let mut gap_size_h = self.gap_size_h[i];
            for j in 1..ncol {
                let b_base = alternate[j - 1];
                let step_diag = last_row[j - 1] + if a_base == b_base { w_match } else { w_mismatch };
                let mut prev_gap = last_row[j] + w_open;
                best_gap_v[j] += w_extend;
                if prev_gap > best_gap_v[j] {
                    best_gap_v[j] = prev_gap;
                    gap_size_v[j] = 1;
                } else {
                    gap_size_v[j] += 1;
                }
                let step_down = best_gap_v[j];
                let kd = gap_size_v[j];
                prev_gap = cur_row[j - 1] + w_open;
                best_gap_h += w_extend;
                if prev_gap > best_gap_h {
                    best_gap_h = prev_gap;
                    gap_size_h = 1;
                } else {
                    gap_size_h += 1;
                }
                let step_right = best_gap_h;
                let ki = gap_size_h;
                if step_diag >= step_down && step_diag >= step_right {
                    cur_row[j] = MATRIX_MIN_CUTOFF.max(step_diag);
                    bt_row[j] = 0;
                } else if step_right >= step_down {
                    cur_row[j] = MATRIX_MIN_CUTOFF.max(step_right);
                    bt_row[j] = -ki;
                } else {
                    cur_row[j] = MATRIX_MIN_CUTOFF.max(step_down);
                    bt_row[j] = kd;
                }
            }
            self.best_gap_h[i] = best_gap_h;
            self.gap_size_h[i] = gap_size_h;
        }
    }

    fn calculate_cigar(&self, ref_len: usize, alt_len: usize, strategy: Overhang) -> SwAlignment {
        let ncol = alt_len + 1;
        let sw = &self.sw;
        let btrack = &self.btrack;
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
                let cur = sw[i * ncol + alt_len];
                if cur >= maxscore {
                    p1 = i as i64;
                    maxscore = cur;
                }
            }
            if strategy != Overhang::LeadingIndel {
                for j in 1..=alt_len {
                    let cur = sw[ref_len * ncol + j];
                    if cur > maxscore
                        || (cur == maxscore && (ref_length - j as i64).abs() < (p1 - p2).abs())
                    {
                        p1 = ref_length;
                        p2 = j as i64;
                        maxscore = cur;
                        segment_length = alt_length - j as i64;
                    }
                }
            }
        }
        let mut lce: Vec<(u8, u32)> = Vec::with_capacity(5);
        if segment_length > 0 && strategy == Overhang::SoftClip {
            lce.push((b'S', segment_length as u32));
            segment_length = 0;
        }
        // states: M=match, I=insertion, D=deletion
        let mut state = b'M';
        loop {
            let btr = btrack[(p1 as usize) * ncol + p2 as usize];
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

/// Column-major variant of the Java Smith-Waterman (SOFTCLIP overhang strategy only) for
/// aligning many alternates against one fixed reference. Column j of the DP matrices depends
/// only on alt[0..j], so columns are reused for the longest common prefix with the previously
/// aligned alternate. Results are bit-identical to `SwAligner::align(.., SoftClip)`.
#[derive(Default)]
pub struct PrefixCachedSw {
    reference: Vec<u8>,
    params: Option<SwParameters>,
    alt: Vec<u8>,
    /// number of alt columns (beyond column 0) currently valid
    valid_cols: usize,
    h: Vec<i32>,
    bgh: Vec<i32>,
    gsh: Vec<i32>,
    bt: Vec<i32>,
}

impl PrefixCachedSw {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_reference(&mut self, reference: &[u8], params: &SwParameters) {
        self.reference.clear();
        self.reference.extend_from_slice(reference);
        self.params = Some(*params);
        self.alt.clear();
        self.valid_cols = 0;
    }

    pub fn align(&mut self, alternate: &[u8]) -> SwAlignment {
        let p = self.params.expect("reference not set");
        assert!(!self.reference.is_empty() && !alternate.is_empty());
        let idx = last_index_of(&self.reference, alternate);
        if idx != -1 {
            return SwAlignment { cigar: vec![(b'M', alternate.len() as u32)], offset: idx };
        }
        let n = self.reference.len();
        let rows = n + 1;
        let m = alternate.len();
        // common prefix with the cached alternate
        let mut lcp = 0usize;
        let lim = self.valid_cols.min(m).min(self.alt.len());
        while lcp < lim && self.alt[lcp] == alternate[lcp] {
            lcp += 1;
        }
        let need = (m + 1) * rows;
        if self.h.len() < need {
            self.h.resize(need, 0);
            self.bgh.resize(need, 0);
            self.gsh.resize(need, 0);
            self.bt.resize(need, 0);
        }
        if lcp == 0 {
            // column 0
            for i in 0..rows {
                self.h[i] = 0;
                self.bgh[i] = LOW_INIT_VALUE;
                self.gsh[i] = 0;
                self.bt[i] = 0;
            }
        }
        let w_open = p.gap_open_penalty;
        let w_extend = p.gap_extend_penalty;
        let w_match = p.match_value;
        let w_mismatch = p.mismatch_penalty;
        let reference = &self.reference[..];
        for j in (lcp + 1)..=m {
            let b_base = alternate[j - 1];
            let (prev_part, cur_part) = self.h.split_at_mut(j * rows);
            let prev_h = &prev_part[(j - 1) * rows..];
            let cur_h = &mut cur_part[..rows];
            let (prev_bgh_part, cur_bgh_part) = self.bgh.split_at_mut(j * rows);
            let prev_bgh = &prev_bgh_part[(j - 1) * rows..];
            let cur_bgh = &mut cur_bgh_part[..rows];
            let (prev_gsh_part, cur_gsh_part) = self.gsh.split_at_mut(j * rows);
            let prev_gsh = &prev_gsh_part[(j - 1) * rows..];
            let cur_gsh = &mut cur_gsh_part[..rows];
            let cur_bt = &mut self.bt[j * rows..(j + 1) * rows];
            cur_h[0] = 0;
            cur_bgh[0] = prev_bgh[0];
            cur_gsh[0] = prev_gsh[0];
            cur_bt[0] = 0;
            let mut best_gap_v = LOW_INIT_VALUE;
            let mut gap_size_v = 0i32;
            let mut above = 0i32; // cur_h[i-1]
            for i in 1..rows {
                let a_base = reference[i - 1];
                let step_diag = prev_h[i - 1] + if a_base == b_base { w_match } else { w_mismatch };
                // vertical gap (step down), from the cell above in this column
                let pg = above + w_open;
                best_gap_v += w_extend;
                if pg > best_gap_v {
                    best_gap_v = pg;
                    gap_size_v = 1;
                } else {
                    gap_size_v += 1;
                }
                let step_down = best_gap_v;
                // horizontal gap (step right), from the cell to the left (previous column)
                let pgh = prev_h[i] + w_open;
                let mut bgh = prev_bgh[i] + w_extend;
                let mut gsh = prev_gsh[i];
                if pgh > bgh {
                    bgh = pgh;
                    gsh = 1;
                } else {
                    gsh += 1;
                }
                cur_bgh[i] = bgh;
                cur_gsh[i] = gsh;
                let step_right = bgh;
                let v;
                let btv;
                if step_diag >= step_down && step_diag >= step_right {
                    v = step_diag;
                    btv = 0;
                } else if step_right >= step_down {
                    v = step_right;
                    btv = -gsh;
                } else {
                    v = step_down;
                    btv = gap_size_v;
                }
                let v = MATRIX_MIN_CUTOFF.max(v);
                cur_h[i] = v;
                cur_bt[i] = btv;
                above = v;
            }
        }
        self.alt.clear();
        self.alt.extend_from_slice(alternate);
        self.valid_cols = m;

        // ---- calculateCigar (SOFTCLIP) on the column-major matrices ----
        let h = &self.h;
        let bt = &self.bt;
        let ref_length = n as i64;
        let alt_length = m as i64;
        let mut p1: i64 = 0;
        let mut p2: i64 = alt_length;
        let mut maxscore = i32::MIN;
        let mut segment_length: i64 = 0;
        for i in 1..=n {
            let cur = h[m * rows + i];
            if cur >= maxscore {
                p1 = i as i64;
                maxscore = cur;
            }
        }
        for j in 1..=m {
            let cur = h[j * rows + n];
            if cur > maxscore || (cur == maxscore && (ref_length - j as i64).abs() < (p1 - p2).abs()) {
                p1 = ref_length;
                p2 = j as i64;
                maxscore = cur;
                segment_length = alt_length - j as i64;
            }
        }
        let mut lce: Vec<(u8, u32)> = Vec::with_capacity(8);
        if segment_length > 0 {
            lce.push((b'S', segment_length as u32));
            segment_length = 0;
        }
        let mut state = b'M';
        loop {
            let btr = bt[(p2 as usize) * rows + p1 as usize];
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
        lce.push((state, segment_length as u32));
        if p2 > 0 {
            lce.push((b'S', p2 as u32));
        }
        lce.reverse();
        SwAlignment { cigar: lce, offset: p1 as i32 }
    }
}

