//! bwa-mem2 index files, memory-mapped (FMI_search.cpp `load_index`, bntseq.cpp).
//!
//! `.bwt.2bit.64`: i64 reference_seq_len (both strands + sentinel), i64 count[5], CP_OCC[(len >> 6) + 1]
//! (per 64 BWT positions: i64 cp_count[4], u64 one_hot_bwt_str[4]), then the suffix array sampled every 8
//! (SA_COMPX 3): i8 sa_ms_byte[(len >> 3) + 1], u32 sa_ls_word[(len >> 3) + 1], i64 sentinel_index.

use std::fs::File;
use std::os::unix::io::AsRawFd;

extern "C" {
    fn mmap(addr: *mut u8, len: usize, prot: i32, flags: i32, fd: i32, off: i64) -> *mut u8;
}
const PROT_READ: i32 = 1;
const MAP_PRIVATE: i32 = 2;

/// A read-only memory map of a whole file (kept for the life of the process).
pub struct Map {
    ptr: *const u8,
    len: usize,
}
unsafe impl Send for Map {}
unsafe impl Sync for Map {}

impl Map {
    pub fn open(path: &str) -> Map {
        let f = File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let len = f.metadata().unwrap().len() as usize;
        let ptr = unsafe { mmap(std::ptr::null_mut(), len, PROT_READ, MAP_PRIVATE, f.as_raw_fd(), 0) };
        assert!(ptr as isize != -1, "{path}: mmap failed");
        Map { ptr, len }
    }
    pub fn bytes(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
    fn i64_at(&self, off: usize) -> i64 {
        i64::from_le_bytes(self.bytes()[off..off + 8].try_into().unwrap())
    }
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct CpOcc {
    pub cp_count: [i64; 4],
    pub one_hot: [u64; 4],
}

pub const CP_SHIFT: i64 = 6;
pub const CP_MASK: i64 = 63;
pub const SA_COMPX: i64 = 3;
pub const SA_COMPX_MASK: i64 = 7;

#[allow(dead_code)] // fields that are never read: kept to mirror bwa-mem2's struct
pub struct Fmi {
    _map: Map,
    pub reference_seq_len: i64,
    /// count[] + 1, as load_index leaves it
    pub count: [i64; 5],
    cp_occ: *const CpOcc,
    n_cp: usize,
    sa_ms_byte: *const i8,
    sa_ls_word: *const u32,
    n_sa: usize,
    pub sentinel_index: i64,
    one_hot_mask: [u64; 64],
}
unsafe impl Send for Fmi {}
unsafe impl Sync for Fmi {}

/// A bi-interval: query [m, n], BWT interval [k, k + s) and its reverse-complement start l.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Smem {
    pub rid: u32,
    pub m: u32,
    pub n: u32,
    pub k: i64,
    pub l: i64,
    pub s: i64,
}

impl Fmi {
    pub fn load(prefix: &str) -> Fmi {
        let map = Map::open(&format!("{prefix}.bwt.2bit.64"));
        let len = map.i64_at(0);
        let mut count = [0i64; 5];
        for (i, c) in count.iter_mut().enumerate() {
            *c = map.i64_at(8 + 8 * i) + 1;
        }
        let n_cp = ((len >> CP_SHIFT) + 1) as usize;
        let cp_off = 48usize;
        let n_sa = ((len >> SA_COMPX) + 1) as usize;
        let ms_off = cp_off + n_cp * 64;
        let ls_off = ms_off + n_sa;
        let sent_off = ls_off + 4 * n_sa;
        assert_eq!(map.bytes().len(), sent_off + 8, "{prefix}.bwt.2bit.64: unexpected size");
        let base = map.bytes().as_ptr();
        let sentinel_index = map.i64_at(sent_off);
        let mut one_hot_mask = [0u64; 64];
        let b = 0x8000_0000_0000_0000u64;
        one_hot_mask[1] = b;
        for i in 2..64 {
            one_hot_mask[i] = (one_hot_mask[i - 1] >> 1) | b;
        }
        unsafe {
            Fmi {
                reference_seq_len: len,
                count,
                cp_occ: base.add(cp_off) as *const CpOcc,
                n_cp,
                sa_ms_byte: base.add(ms_off) as *const i8,
                sa_ls_word: base.add(ls_off) as *const u32,
                n_sa,
                sentinel_index,
                one_hot_mask,
                _map: map,
            }
        }
    }

    /// Raw parts for the GPU: checkpoints, SA sample bytes, SA sample words, their counts
    #[cfg_attr(not(gpu), allow(dead_code))] // GPU build only
    pub fn raw(&self) -> (*const CpOcc, usize, *const i8, *const u32, usize) {
        (self.cp_occ, self.n_cp, self.sa_ms_byte, self.sa_ls_word, self.n_sa)
    }

    #[inline]
    fn cp(&self, i: i64) -> &CpOcc {
        debug_assert!((i as usize) < self.n_cp);
        unsafe { &*self.cp_occ.add(i as usize) }
    }

    /// GET_OCC: occurrences of `c` in BWT[0, pp)
    #[inline]
    fn occ(&self, pp: i64, c: usize) -> i64 {
        let cp = self.cp(pp >> CP_SHIFT);
        let y = (pp & CP_MASK) as usize;
        cp.cp_count[c] + (cp.one_hot[c] & self.one_hot_mask[y]).count_ones() as i64
    }

    /// `FMI_search::backwardExt`
    #[inline]
    pub fn backward_ext(&self, mut smem: Smem, a: u8) -> Smem {
        let mut k = [0i64; 4];
        let mut s = [0i64; 4];
        let sp = smem.k;
        let ep = smem.k + smem.s;
        for b in 0..4 {
            let occ_sp = self.occ(sp, b);
            let occ_ep = self.occ(ep, b);
            k[b] = self.count[b] + occ_sp;
            s[b] = occ_ep - occ_sp;
        }
        let sentinel_offset = if smem.k <= self.sentinel_index && smem.k + smem.s > self.sentinel_index { 1 } else { 0 };
        let mut l = [0i64; 4];
        l[3] = smem.l + sentinel_offset;
        l[2] = l[3] + s[3];
        l[1] = l[2] + s[2];
        l[0] = l[1] + s[1];
        let a = a as usize;
        smem.k = k[a];
        smem.l = l[a];
        smem.s = s[a];
        smem
    }

    #[inline]
    fn sa_sample(&self, i: i64) -> i64 {
        debug_assert!((i as usize) < self.n_sa);
        unsafe { ((*self.sa_ms_byte.add(i as usize) as i64) << 32) + *self.sa_ls_word.add(i as usize) as i64 }
    }

    /// The suffix-array entry at `pos`, as `get_sa_entries_prefetch` / `call_one_step` compute it: LF-walk back to a
    /// sampled position (every 8th), adding the steps. Reaching the sentinel ('$', no base in the one-hot rows)
    /// gives 0, not the steps walked (call_one_step's `sa_entry = 0`, unlike get_sa_entry_compressed).
    pub fn sa_entry(&self, pos: i64) -> i64 {
        if pos & SA_COMPX_MASK == 0 {
            return self.sa_sample(pos >> SA_COMPX);
        }
        let mut sp = pos;
        let mut offset = 0i64;
        loop {
            let cp = self.cp(sp >> CP_SHIFT);
            let y = CP_MASK - (sp & CP_MASK);
            let b = if (cp.one_hot[0] >> y) & 1 != 0 {
                0
            } else if (cp.one_hot[1] >> y) & 1 != 0 {
                1
            } else if (cp.one_hot[2] >> y) & 1 != 0 {
                2
            } else if (cp.one_hot[3] >> y) & 1 != 0 {
                3
            } else {
                return 0;
            };
            sp = self.count[b] + self.occ(sp, b);
            offset += 1;
            if sp & SA_COMPX_MASK == 0 {
                return self.sa_sample(sp >> SA_COMPX) + offset;
            }
        }
    }
}

/// One reference sequence of `.ann` (bntann1_t).
#[allow(dead_code)] // fields that are never read: kept to mirror bwa-mem2's struct
#[derive(Clone, Debug)]
pub struct Ann {
    pub offset: i64,
    pub len: i32,
    pub n_ambs: i32,
    pub gi: u32,
    pub is_alt: bool,
    pub name: String,
    pub anno: String,
}

/// `.amb` hole (bntamb1_t)
#[allow(dead_code)] // fields that are never read: kept to mirror bwa-mem2's struct
#[derive(Clone, Copy, Debug)]
pub struct Amb {
    pub offset: i64,
    pub len: i32,
    pub amb: u8,
}

/// bntseq_t: reference annotation plus the 2-bit packed forward strand (`.pac`).
#[allow(dead_code)] // fields that are never read: kept to mirror bwa-mem2's struct
pub struct Bns {
    pub l_pac: i64,
    pub seed: u32,
    pub anns: Vec<Ann>,
    pub ambs: Vec<Amb>,
    pac: Map,
    /// anns[i].offset, packed (pos2rid's binary search touched one 80-byte Ann per step)
    offs: Vec<i64>,
    /// per 2^RID_SHIFT-base block of the forward strand: the reference holding all of it, or -1 when a reference
    /// starts inside it (pos2rid searches `offs` then)
    rid_tab: Vec<i32>,
}

const RID_SHIFT: u32 = 16;

impl Bns {
    pub fn load(prefix: &str) -> Bns {
        let ann = std::fs::read_to_string(format!("{prefix}.ann")).unwrap_or_else(|e| panic!("{prefix}.ann: {e}"));
        let mut lines = ann.lines();
        let h: Vec<i64> = lines.next().unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect();
        let (l_pac, n_seqs, seed) = (h[0], h[1] as usize, h[2] as u32);
        let mut anns = Vec::with_capacity(n_seqs);
        for _ in 0..n_seqs {
            // "gi name anno" then "offset len n_ambs"
            let l1 = lines.next().unwrap();
            let mut it = l1.splitn(3, ' ');
            let gi: u32 = it.next().unwrap().parse().unwrap();
            let name = it.next().unwrap().to_string();
            // bns_restore_core: the rest of the line after one space; "(null)" or nothing -> ""
            let anno = it.next().filter(|a| !a.is_empty() && *a != "(null)").unwrap_or("").to_string();
            let f: Vec<i64> = lines.next().unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect();
            anns.push(Ann { offset: f[0], len: f[1] as i32, n_ambs: f[2] as i32, gi, is_alt: false, name, anno });
        }
        let amb = std::fs::read_to_string(format!("{prefix}.amb")).unwrap_or_else(|e| panic!("{prefix}.amb: {e}"));
        let mut lines = amb.lines();
        let h: Vec<i64> = lines.next().unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect();
        assert_eq!(h[0], l_pac, ".amb and .ann disagree on l_pac");
        let mut ambs = Vec::with_capacity(h[2] as usize);
        for _ in 0..h[2] {
            let f: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
            ambs.push(Amb { offset: f[0].parse().unwrap(), len: f[1].parse().unwrap(), amb: f[2].as_bytes()[0] });
        }
        // .alt: names of ALT contigs (bwa_idx_load_from_disk -> bns_restore with the .alt file)
        // bns_restore reads it by character: a name counts when a tab, newline or CR ends it (the rest of that line is
        // skipped; '@' lines are headers); a name at the end of the file without one is not looked up
        if let Ok(alt) = std::fs::read(format!("{prefix}.alt")) {
            let names: std::collections::HashMap<&[u8], usize> = anns.iter().enumerate().map(|(i, a)| (a.name.as_bytes(), i)).collect();
            let mut alt_ids = Vec::new();
            let mut i = 0;
            let mut st = 0;
            while i < alt.len() {
                let c = alt[i];
                if c == b'\t' || c == b'\n' || c == b'\r' {
                    let tok = &alt[st..i];
                    if tok.first() != Some(&b'@') {
                        if let Some(&k) = names.get(tok) {
                            alt_ids.push(k);
                        }
                    }
                    while i < alt.len() && alt[i] != b'\n' {
                        i += 1;
                    }
                    i += 1;
                    st = i;
                } else {
                    i += 1;
                }
            }
            for k in alt_ids {
                anns[k].is_alt = true;
            }
        }
        let pac = Map::open(&format!("{prefix}.pac"));
        let offs = anns.iter().map(|a| a.offset).collect();
        let mut b = Bns { l_pac, seed, anns, ambs, pac, offs, rid_tab: Vec::new() };
        // pos2rid is the last reference starting at or before pos: non-decreasing in pos, so a block whose first and
        // last positions give the same reference gives it everywhere
        b.rid_tab = (0..(l_pac + (1 << RID_SHIFT) - 1) >> RID_SHIFT)
            .map(|k| {
                let (s, e) = (k << RID_SHIFT, (((k + 1) << RID_SHIFT) - 1).min(l_pac - 1));
                let r = b.pos2rid_search(s);
                if b.pos2rid_search(e) == r { r } else { -1 }
            })
            .collect();
        b
    }

    /// The 2-bit packed forward strand (.pac)
    pub fn pac_bytes(&self) -> &[u8] {
        self.pac.bytes()
    }

    /// Base `i` (0..l_pac) of the forward strand: _get_pac(pac, i)
    #[cfg_attr(not(gpu), allow(dead_code))] // GPU build only
    #[inline]
    pub fn pac_base(&self, i: i64) -> u8 {
        (self.pac.bytes()[(i >> 2) as usize] >> ((!(i as u64) & 3) << 1)) & 3
    }

    /// `bns_pos2rid`
    #[inline]
    pub fn pos2rid(&self, pos_f: i64) -> i32 {
        if pos_f >= self.l_pac {
            return -1;
        }
        if pos_f >= 0 {
            let r = self.rid_tab[(pos_f >> RID_SHIFT) as usize];
            if r >= 0 {
                return r;
            }
        }
        self.pos2rid_search(pos_f)
    }

    /// `bns_pos2rid`'s binary search
    fn pos2rid_search(&self, pos_f: i64) -> i32 {
        if pos_f >= self.l_pac {
            return -1;
        }
        let offs = &self.offs;
        let (mut left, mut right) = (0usize, offs.len());
        let mut mid = 0usize;
        while left < right {
            mid = (left + right) >> 1;
            if pos_f >= offs[mid] {
                if mid == offs.len() - 1 {
                    break;
                }
                if pos_f < offs[mid + 1] {
                    break;
                }
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        mid as i32
    }

    /// `bns_intv2rid`: the reference sequence holding [rb, re), -1 if it crosses sequences or strands
    pub fn intv2rid(&self, rb: i64, re: i64) -> i32 {
        if rb < self.l_pac && re > self.l_pac {
            return -2;
        }
        assert!(rb <= re);
        let rid_b = self.pos2rid(self.depos(rb).0);
        let rid_e = if rb < re { self.pos2rid(self.depos(re - 1).0) } else { rid_b };
        if rid_b == rid_e {
            rid_b
        } else {
            -1
        }
    }

    /// `bns_depos`: forward-strand position and whether `pos` was on the reverse strand
    #[inline]
    pub fn depos(&self, pos: i64) -> (i64, bool) {
        if pos >= self.l_pac {
            ((self.l_pac << 1) - 1 - pos, true)
        } else {
            (pos, false)
        }
    }
}
