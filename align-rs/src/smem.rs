//! Seeding (bwamem.cpp `mem_collect_smem`, FMI_search.cpp `getSMEMsOnePosOneThread`,
//! `getSMEMsAllPosOneThread`, `bwtSeedStrategyAllPosOneThread`) for one read: SMEMs from every position, re-seeding
//! inside long SMEMs with few hits, then the LAST-like third round; sorted by (m, n) as `ks_introsort(mem_intv1)`
//! leaves them (equal (m, n) are the same interval, so the sort's tie order does not matter).

use crate::index::{Fmi, Smem};

pub struct Opt {
    pub min_seed_len: i32,
    pub split_factor: f64,
    pub split_width: i64,
    pub max_mem_intv: i64,
}

/// getSMEMsOnePosOneThread for one read and one start `x`; returns next_x.
fn one_pos(fmi: &Fmi, q: &[u8], rid: u32, x: usize, min_intv: i64, min_seed_len: i32, out: &mut Vec<Smem>, prev: &mut Vec<Smem>) -> usize {
    let len = q.len();
    let mut next_x = x + 1;
    let mut a = q[x];
    if a >= 4 {
        return next_x;
    }
    let c = &fmi.count;
    let mut smem = Smem { rid, m: x as u32, n: x as u32, k: c[a as usize], l: c[3 - a as usize], s: c[a as usize + 1] - c[a as usize] };
    prev.clear();
    let mut j = x + 1;
    while j < len {
        a = q[j];
        next_x = j + 1;
        if a < 4 {
            // forward extension = backward extension on the reverse complement
            let mut sm_ = smem;
            sm_.k = smem.l;
            sm_.l = smem.k;
            let ns_ = fmi.backward_ext(sm_, 3 - a);
            let mut ns = ns_;
            ns.k = ns_.l;
            ns.l = ns_.k;
            ns.n = j as u32;
            if ns.s != smem.s {
                prev.push(smem);
            }
            if ns.s < min_intv {
                next_x = j;
                break;
            }
            smem = ns;
        } else {
            break;
        }
        j += 1;
    }
    if smem.s >= min_intv {
        prev.push(smem);
    }
    prev.reverse();
    // backward search
    let mut jj = x as i64 - 1;
    while jj >= 0 {
        let j = jj as usize;
        a = q[j];
        if a > 3 {
            break;
        }
        let mut num_curr = 0usize;
        // `int curr_s` in bwa-mem2: the interval size as a 32-bit int
        let mut curr_s: i64 = -1;
        let mut p = 0usize;
        while p < prev.len() {
            let sm = prev[p];
            let mut ns = fmi.backward_ext(sm, a);
            ns.m = j as u32;
            if ns.s < min_intv && (sm.n - sm.m + 1) as i32 >= min_seed_len {
                out.push(sm);
                break;
            }
            if ns.s >= min_intv && ns.s != curr_s {
                curr_s = ns.s as i32 as i64;
                prev[num_curr] = ns;
                num_curr += 1;
                break;
            }
            p += 1;
        }
        p += 1;
        while p < prev.len() {
            let sm = prev[p];
            let mut ns = fmi.backward_ext(sm, a);
            ns.m = j as u32;
            if ns.s >= min_intv && ns.s != curr_s {
                curr_s = ns.s as i32 as i64;
                prev[num_curr] = ns;
                num_curr += 1;
            }
            p += 1;
        }
        prev.truncate(num_curr);
        if num_curr == 0 {
            break;
        }
        jj -= 1;
    }
    if !prev.is_empty() {
        let sm = prev[0];
        if (sm.n - sm.m + 1) as i32 >= min_seed_len {
            out.push(sm);
        }
        prev.clear();
    }
    next_x
}

/// bwtSeedStrategyAllPosOneThread for one read
fn seed_strategy(fmi: &Fmi, q: &[u8], rid: u32, max_intv: i64, min_seed_len: i32, out: &mut Vec<Smem>) {
    let len = q.len();
    let c = &fmi.count;
    let mut x = 0usize;
    while x < len {
        let mut next_x = x + 1;
        let mut a = q[x];
        if a < 4 {
            let mut smem = Smem { rid, m: x as u32, n: x as u32, k: c[a as usize], l: c[3 - a as usize], s: c[a as usize + 1] - c[a as usize] };
            let mut j = x + 1;
            while j < len {
                next_x = j + 1;
                a = q[j];
                if a < 4 {
                    let mut sm_ = smem;
                    sm_.k = smem.l;
                    sm_.l = smem.k;
                    let ns_ = fmi.backward_ext(sm_, 3 - a);
                    let mut ns = ns_;
                    ns.k = ns_.l;
                    ns.l = ns_.k;
                    ns.n = j as u32;
                    smem = ns;
                    if smem.s < max_intv && (smem.n - smem.m + 1) as i32 >= min_seed_len {
                        if smem.s > 0 {
                            out.push(smem);
                        }
                        break;
                    }
                } else {
                    break;
                }
                j += 1;
            }
        }
        x = next_x;
    }
}

/// All seeds of one read (encoded 0..4), sorted by (m, n).
pub fn collect(fmi: &Fmi, opt: &Opt, q: &[u8], rid: u32) -> Vec<Smem> {
    let mut out = Vec::new();
    let mut prev = Vec::with_capacity(q.len() + 1);
    // round 1: SMEMs from every position (min interval size 1)
    let mut x = 0usize;
    while x < q.len() {
        x = one_pos(fmi, q, rid, x, 1, opt.min_seed_len, &mut out, &mut prev);
    }
    // round 2: re-seed from the middle of long SMEMs with few occurrences
    let split_len = (opt.min_seed_len as f64 * opt.split_factor + 0.499) as i32;
    let n1 = out.len();
    for i in 0..n1 {
        let p = out[i];
        let (start, end) = (p.m as i32, p.n as i32 + 1);
        if end - start < split_len || p.s > opt.split_width {
            continue;
        }
        let pos = ((end + start) >> 1) as usize;
        let min_intv = (p.s + 1) as i32 as i64;
        one_pos(fmi, q, rid, pos, min_intv, opt.min_seed_len, &mut out, &mut prev);
    }
    // round 3
    if opt.max_mem_intv > 0 {
        seed_strategy(fmi, q, rid, opt.max_mem_intv, opt.min_seed_len + 1, &mut out);
    }
    out.sort_by_key(|s| ((s.m as u64) << 32) | s.n as u64);
    out
}
