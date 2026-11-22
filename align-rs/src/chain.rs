//! Chaining (bwamem.cpp `mem_chain_seeds`, `test_and_merge`, `mem_chain_flt`, `mem_flt_chained_seeds`) for one read.
//! Chains are kept in klib's B-tree keyed by position (kbtree.h, ported in `KbTree`): with equal positions, where a
//! new chain goes and which chain a seed is offered to depend on the tree's shape, so the tree itself is ported.

use crate::index::{Bns, Smem};

#[derive(Clone, Copy, Debug, Default)]
pub struct Seed {
    pub rbeg: i64,
    pub qbeg: i32,
    pub len: i32,
    pub score: i32,
}

#[allow(dead_code)] // fields that are never read: kept to mirror bwa-mem2's struct
#[derive(Clone, Debug, Default)]
pub struct Chain {
    pub seqid: i32,
    pub n_alloc: i32,
    pub first: i32,
    pub rid: i32,
    pub w: i32,
    pub kept: i32,
    pub is_alt: bool,
    pub frac_rep: f32,
    pub pos: i64,
    pub seeds: Vec<Seed>,
}

pub struct ChainOpt {
    pub w: i64,
    pub max_chain_gap: i64,
    pub max_occ: i64,
    pub mask_level: f32,
    pub drop_ratio: f32,
    pub min_seed_len: i32,
    pub min_chain_weight: i32,
    pub max_chain_extend: i32,
    pub a: i32,
}

// ------------------------------------------------------------------------------------------------
// kbtree.h, KBTREE_INIT(chn, mem_chain_t, chain_cmp) with kb_init(chn, KB_DEFAULT_SIZE + 8):
// sizeof(mem_chain_t) = 48 -> t = ((520 - 4 - 8) / (8 + 48) + 1) >> 1 = 5, at most 2t - 1 = 9 keys per node.
// Keys are indexes into the chain arena; the comparison is by chain position.
// ------------------------------------------------------------------------------------------------

const KB_T: usize = 5;
const KB_N: usize = 2 * KB_T - 1;

/// A node: its keys (chain indexes) with each key's chain position beside it (chains never move), and children
#[derive(Clone, Copy)]
struct Node {
    internal: bool,
    n: usize,
    keys: [u32; KB_N],
    pos: [i64; KB_N],
    ptrs: [u32; KB_N + 1],
}

impl Node {
    fn new(internal: bool) -> Node {
        Node { internal, n: 0, keys: [0; KB_N], pos: [0; KB_N], ptrs: [0; KB_N + 1] }
    }
    fn insert_key(&mut self, i: usize, key: u32, pos: i64) {
        self.keys.copy_within(i..self.n, i + 1);
        self.pos.copy_within(i..self.n, i + 1);
        self.keys[i] = key;
        self.pos[i] = pos;
        self.n += 1;
    }
}

pub struct KbTree {
    nodes: Vec<Node>,
    root: usize,
    pub n_keys: usize,
}

#[inline]
fn cmp(a: i64, b: i64) -> i32 {
    ((b < a) as i32) - ((a < b) as i32)
}

impl KbTree {
    pub fn new() -> KbTree {
        KbTree { nodes: vec![Node::new(false)], root: 0, n_keys: 0 }
    }

    /// __kb_getp_aux: (index, r)
    #[inline]
    fn getp_aux(&self, x: usize, k: i64) -> (i32, i32) {
        let node = &self.nodes[x];
        let n = node.n;
        if n == 0 {
            return (-1, 0);
        }
        let (mut begin, mut end) = (0usize, n);
        while begin < end {
            let mid = (begin + end) >> 1;
            if cmp(node.pos[mid], k) < 0 {
                begin = mid + 1;
            } else {
                end = mid;
            }
        }
        if begin == n {
            return (n as i32 - 1, 1);
        }
        let r = cmp(k, node.pos[begin]);
        let i = if r < 0 { begin as i32 - 1 } else { begin as i32 };
        (i, r)
    }

    /// kb_intervalp: the chain keys lower (<= k) and upper (>= k)
    pub fn interval(&self, k: i64) -> (Option<usize>, Option<usize>) {
        let (mut lower, mut upper) = (None, None);
        let mut x = self.root;
        loop {
            let (i, r) = self.getp_aux(x, k);
            let node = &self.nodes[x];
            if i >= 0 && r == 0 {
                return (Some(node.keys[i as usize] as usize), Some(node.keys[i as usize] as usize));
            }
            if i >= 0 {
                lower = Some(node.keys[i as usize] as usize);
            }
            if (i as i64) < node.n as i64 - 1 {
                upper = Some(node.keys[(i + 1) as usize] as usize);
            }
            if !node.internal {
                return (lower, upper);
            }
            x = node.ptrs[(i + 1) as usize] as usize;
        }
    }

    /// __kb_split: split child y = x.ptrs[i] (full) around its middle key
    fn split(&mut self, x: usize, i: usize, y: usize) {
        let yn = self.nodes[y];
        let mut z = Node::new(yn.internal);
        z.n = KB_N - KB_T;
        z.keys[..z.n].copy_from_slice(&yn.keys[KB_T..KB_N]);
        z.pos[..z.n].copy_from_slice(&yn.pos[KB_T..KB_N]);
        if yn.internal {
            z.ptrs[..KB_T].copy_from_slice(&yn.ptrs[KB_T..KB_N + 1]);
        }
        let (mid, mid_pos) = (yn.keys[KB_T - 1], yn.pos[KB_T - 1]);
        self.nodes[y].n = KB_T - 1;
        let zi = self.nodes.len() as u32;
        self.nodes.push(z);
        let xn = &mut self.nodes[x];
        let np = xn.n + 1;
        xn.ptrs.copy_within(i + 1..np, i + 2);
        xn.ptrs[i + 1] = zi;
        xn.insert_key(i, mid, mid_pos);
    }

    fn putp_aux(&mut self, x: usize, key: usize, k: i64) {
        if !self.nodes[x].internal {
            let (i, _) = self.getp_aux(x, k);
            self.nodes[x].insert_key((i + 1) as usize, key as u32, k);
        } else {
            let mut i = (self.getp_aux(x, k).0 + 1) as usize;
            let child = self.nodes[x].ptrs[i] as usize;
            if self.nodes[child].n == KB_N {
                self.split(x, i, child);
                if cmp(k, self.nodes[x].pos[i]) > 0 {
                    i += 1;
                }
            }
            let child = self.nodes[x].ptrs[i] as usize;
            self.putp_aux(child, key, k);
        }
    }

    /// kb_putp: `key` indexes a chain at position `pos`
    pub fn put(&mut self, key: usize, pos: i64) {
        self.n_keys += 1;
        let r = self.root;
        if self.nodes[r].n == KB_N {
            let s = self.nodes.len();
            let mut root = Node::new(true);
            root.ptrs[0] = r as u32;
            self.nodes.push(root);
            self.root = s;
            self.split(s, 0, r);
        }
        let root = self.root;
        self.putp_aux(root, key, pos);
    }

    /// __kb_traverse: keys in order
    pub fn traverse(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.n_keys);
        self.walk(self.root, &mut out);
        out
    }
    fn walk(&self, x: usize, out: &mut Vec<usize>) {
        let node = &self.nodes[x];
        for i in 0..=node.n {
            if node.internal {
                self.walk(node.ptrs[i] as usize, out);
            }
            if i < node.n {
                out.push(node.keys[i] as usize);
            }
        }
    }
}

const SEEDS_PER_CHAIN: i32 = 1; // only the allocation size bwa-mem2 starts a chain with (not used for output)

/// test_and_merge: 1 if the seed joins (or is contained in) chain `c`
fn test_and_merge(opt: &ChainOpt, l_pac: i64, c: &mut Chain, p: &Seed, seed_rid: i32) -> bool {
    let last = *c.seeds.last().unwrap();
    let qend = last.qbeg as i64 + last.len as i64;
    let rend = last.rbeg + last.len as i64;
    if seed_rid != c.rid {
        return false;
    }
    let s0 = c.seeds[0];
    if p.qbeg >= s0.qbeg && p.qbeg as i64 + p.len as i64 <= qend && p.rbeg >= s0.rbeg && p.rbeg + p.len as i64 <= rend {
        return true; // contained seed
    }
    if (last.rbeg < l_pac || s0.rbeg < l_pac) && p.rbeg >= l_pac {
        return false; // different strand
    }
    let x = p.qbeg as i64 - last.qbeg as i64;
    let y = p.rbeg - last.rbeg;
    if y >= 0 && x - y <= opt.w && y - x <= opt.w && x - (last.len as i64) < opt.max_chain_gap && y - (last.len as i64) < opt.max_chain_gap {
        c.seeds.push(*p);
        return true;
    }
    false
}

/// mem_chain_seeds for one read (its sorted seeds `smems`); `l_seq` its length, `seqid` its index in the chunk.
/// The suffix-array positions mem_chain_seeds looks up for `smems`, in its order (each seed's interval, stepped
/// so at most max_occ occurrences).
pub fn sa_positions(opt: &ChainOpt, smems: &[Smem], out: &mut Vec<i64>) {
    for p in smems {
        let step = if p.s > opt.max_occ { p.s / opt.max_occ } else { 1 };
        let (mut k, mut count) = (0i64, 0i64);
        while k < p.s && count < opt.max_occ {
            out.push(p.k + k);
            k += step;
            count += 1;
        }
    }
}

/// The number of positions sa_positions gives for `smems`
pub fn sa_count(opt: &ChainOpt, smems: &[Smem]) -> usize {
    let mut n = 0usize;
    for p in smems {
        let step = if p.s > opt.max_occ { p.s / opt.max_occ } else { 1 };
        // k = 0, step, 2 step, ... < s, at most max_occ of them
        let c = (p.s + step - 1) / step;
        n += c.min(opt.max_occ) as usize;
    }
    n
}

/// `coords`: the SA entries of sa_positions(smems), in that order
pub fn chain_seeds(bns: &Bns, opt: &ChainOpt, smems: &[Smem], coords: &[i64], l_seq: i32, seqid: i32) -> Vec<Chain> {
    let l_pac = bns.l_pac;
    // frac_rep
    let (mut b, mut e, mut l_rep) = (0i32, 0i32, 0i32);
    for p in smems {
        let (sb, se) = (p.m as i32, p.n as i32 + 1);
        if p.s <= opt.max_occ {
            continue;
        }
        if sb > e {
            l_rep += e - b;
            b = sb;
            e = se;
        } else {
            e = e.max(se);
        }
    }
    l_rep += e - b;
    let mut chains: Vec<Chain> = Vec::with_capacity(coords.len().min(1024));
    let mut tree = KbTree::new();
    let mut ci = 0usize;
    for p in smems {
        let slen = (p.n + 1 - p.m) as i32;
        let step = if p.s > opt.max_occ { p.s / opt.max_occ } else { 1 };
        let (mut k, mut count) = (0i64, 0i64);
        while k < p.s && count < opt.max_occ {
            let rbeg = coords[ci];
            ci += 1;
            k += step;
            count += 1;
            let s = Seed { rbeg, qbeg: p.m as i32, len: slen, score: slen };
            let rid = bns.intv2rid(s.rbeg, s.rbeg + s.len as i64);
            if rid < 0 {
                continue;
            }
            let to_add = if tree.n_keys > 0 {
                let (lower, _) = tree.interval(s.rbeg);
                match lower {
                    None => true,
                    Some(ci) => !test_and_merge(opt, l_pac, &mut chains[ci], &s, rid),
                }
            } else {
                true
            };
            if to_add {
                chains.push(Chain {
                    seqid,
                    n_alloc: SEEDS_PER_CHAIN,
                    first: 0,
                    rid,
                    w: 0,
                    kept: 0,
                    is_alt: bns.anns[rid as usize].is_alt,
                    frac_rep: 0.0,
                    pos: s.rbeg,
                    seeds: vec![s],
                });
                let key = chains.len() - 1;
                tree.put(key, s.rbeg);
            }
        }
    }
    let frac = l_rep as f32 / l_seq as f32;
    // in the tree's order (each chain once: the moved-out slot is left empty)
    tree.traverse()
        .into_iter()
        .map(|i| {
            let mut c = std::mem::take(&mut chains[i]);
            c.frac_rep = frac;
            c
        })
        .collect()
}

/// mem_chain_weight
pub fn chain_weight(c: &Chain) -> i32 {
    let (mut end, mut w) = (0i64, 0i32);
    for s in &c.seeds {
        let (qb, qe) = (s.qbeg as i64, s.qbeg as i64 + s.len as i64);
        if qb >= end {
            w += s.len;
        } else if qe > end {
            w += (qe - end) as i32;
        }
        end = end.max(qe);
    }
    let tmp = w;
    let (mut end, mut w) = (0i64, 0i32);
    for s in &c.seeds {
        let (rb, re) = (s.rbeg, s.rbeg + s.len as i64);
        if rb >= end {
            w += s.len;
        } else if re > end {
            w += (re - end) as i32;
        }
        end = end.max(re);
    }
    let w = w.min(tmp);
    if w < 1 << 30 {
        w
    } else {
        (1 << 30) - 1
    }
}

fn chn_beg(c: &Chain) -> i32 {
    c.seeds[0].qbeg
}
fn chn_end(c: &Chain) -> i32 {
    let l = c.seeds.last().unwrap();
    l.qbeg + l.len
}

/// mem_chain_flt for one read's chains (one seqid, so one range)
pub fn chain_flt(opt: &ChainOpt, mut a: Vec<Chain>) -> Vec<Chain> {
    if a.is_empty() {
        return a;
    }
    for c in a.iter_mut() {
        c.first = -1;
        c.kept = 0;
        c.w = chain_weight(c);
    }
    a.retain(|c| c.w >= opt.min_chain_weight);
    let n = a.len();
    if n == 0 {
        return a;
    }
    // ks_introsort(mem_flt): by weight, descending, on indexes (the chains move as wholes)
    let mut idx: Vec<usize> = (0..n).collect();
    {
        let w: Vec<i32> = a.iter().map(|c| c.w).collect();
        crate::ksort::introsort(&mut idx, |x, y| w[*x] > w[*y]);
    }
    let mut taken: Vec<Option<Chain>> = a.into_iter().map(Some).collect();
    let mut a: Vec<Chain> = idx.into_iter().map(|i| taken[i].take().unwrap()).collect();
    // the pairwise test below on flat per-chain arrays (a repetitive read has thousands of chains, and reaching a
    // chain's ends through its seed list missed the cache for every pair)
    let beg: Vec<i32> = a.iter().map(chn_beg).collect();
    let end: Vec<i32> = a.iter().map(chn_end).collect();
    let wt: Vec<i32> = a.iter().map(|c| c.w).collect();
    let alt: Vec<bool> = a.iter().map(|c| c.is_alt).collect();
    let mut first: Vec<i32> = vec![-1; n];
    let mut kept: Vec<i32> = vec![0; n];
    kept[0] = 3;
    let mut kept_idx: Vec<usize> = vec![0];
    for i in 1..n {
        let mut large_ovlp = false;
        let mut k = 0;
        let (bi, ei, wi, alti) = (beg[i], end[i], wt[i], alt[i]);
        let li = ei - bi;
        while k < kept_idx.len() {
            let j = kept_idx[k];
            let b_max = beg[j].max(bi);
            let e_min = end[j].min(ei);
            if e_min > b_max && (!alt[j] || alti) {
                let lj = end[j] - beg[j];
                let min_l = li.min(lj);
                if ((e_min - b_max) as f32) >= min_l as f32 * opt.mask_level && min_l < opt.max_chain_gap as i32 {
                    large_ovlp = true;
                    if first[j] < 0 {
                        first[j] = i as i32;
                    }
                    if (wi as f32) < wt[j] as f32 * opt.drop_ratio && wt[j] - wi >= opt.min_seed_len << 1 {
                        break;
                    }
                }
            }
            k += 1;
        }
        if k == kept_idx.len() {
            kept_idx.push(i);
            kept[i] = if large_ovlp { 2 } else { 3 };
        }
    }
    for &ci in &kept_idx {
        let f = first[ci];
        if f >= 0 {
            kept[f as usize] = 1;
        }
    }
    for (c, (&f, &k)) in a.iter_mut().zip(first.iter().zip(&kept)) {
        c.first = f;
        c.kept = k;
    }
    let (mut i, mut k) = (0usize, 0i32);
    while i < n {
        if a[i].kept == 0 || a[i].kept == 3 {
            i += 1;
            continue;
        }
        k += 1;
        if k >= opt.max_chain_extend {
            break;
        }
        i += 1;
    }
    while i < n {
        if a[i].kept < 3 {
            a[i].kept = 0;
        }
        i += 1;
    }
    a.retain(|c| c.kept != 0);
    a
}

const MEM_SHORT_EXT: i64 = 50;
const MEM_SHORT_LEN: i64 = 200;

/// mem_flt_chained_seeds for one read's chains (`query`: its bases, 0..4). Reads where min_l > 0.05 * l_query
/// (with the default -W 0: shorter than 725 bp) are left alone; for longer ones each seed is re-scored by a local
/// alignment of its surroundings (mem_seed_sw), seeds scoring under min_HSP_score are dropped and the others take
/// that score (extension orders seeds by it).
pub fn flt_chained_seeds(opt: &ChainOpt, e: &crate::ext::ExtOpt, bns: &Bns, query: &[u8], chains: &mut [Chain]) {
    let l_query = query.len() as i32;
    let min_l: f64 = if opt.min_chain_weight != 0 { 1.1f32 as f64 * opt.min_chain_weight as f64 } else { 5.5f32 as f64 * (l_query as f64).ln() };
    if min_l > (0.05f32 * l_query as f32) as f64 {
        return;
    }
    let min_hsp_score = (opt.a as f64 * min_l + 0.499) as i32;
    for c in chains.iter_mut() {
        c.seeds.retain_mut(|s| {
            let sc = seed_sw(e, bns, query, s);
            if sc < 0 || sc >= min_hsp_score {
                s.score = if sc < 0 { s.len * opt.a } else { sc };
                true
            } else {
                false
            }
        });
    }
}

/// mem_seed_sw: the local-alignment score of a short seed's surroundings (50 bp each side), or -1 when the seed
/// or its window is long enough not to need it
fn seed_sw(e: &crate::ext::ExtOpt, bns: &Bns, query: &[u8], s: &Seed) -> i32 {
    let l_pac = bns.l_pac;
    if s.len as i64 >= MEM_SHORT_LEN {
        return -1;
    }
    let (mut qb, mut qe) = (s.qbeg as i64, s.qbeg as i64 + s.len as i64);
    let (mut rb, mut re) = (s.rbeg, s.rbeg + s.len as i64);
    let mid = (rb + re) >> 1;
    qb = (qb - MEM_SHORT_EXT).max(0);
    qe = (qe + MEM_SHORT_EXT).min(query.len() as i64);
    rb = (rb - MEM_SHORT_EXT).max(0);
    re = (re + MEM_SHORT_EXT).min(l_pac << 1);
    if rb < l_pac && l_pac < re {
        if mid < l_pac {
            re = l_pac;
        } else {
            rb = l_pac;
        }
    }
    if qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN {
        return -1;
    }
    // bns_fetch_seq: the window clipped to mid's contig (on mid's strand)
    let rid = bns.pos2rid(bns.depos(mid).0);
    let (rb, re, _) = crate::ext::fetch_window(bns, rb, mid, re, rid);
    let rseq = crate::ext::bns_get_seq(bns, rb, re);
    let x = crate::ksw::ksw_align2(&query[qb as usize..qe as usize], &rseq, &e.mat, e.o_del, e.e_del, e.o_ins, e.e_ins, crate::ksw::KSW_XSTART);
    x.score
}
