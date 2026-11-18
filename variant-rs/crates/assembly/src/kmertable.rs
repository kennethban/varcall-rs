//! Small open-addressing hash table keyed by k-mers stored in a shared byte arena
//! (key = offset into the arena, fixed k). Linear probing with backward-shift deletion.

#[inline]
pub fn hash_kmer(s: &[u8]) -> u64 {
    // FxHash-style over 8-byte words
    const K: u64 = 0x517c_c1b7_2722_0a95;
    let mut h: u64 = s.len() as u64;
    let mut chunks = s.chunks_exact(8);
    for c in &mut chunks {
        let w = u64::from_le_bytes(c.try_into().unwrap());
        h = (h.rotate_left(5) ^ w).wrapping_mul(K);
    }
    let rem = chunks.remainder();
    if !rem.is_empty() {
        let mut buf = [0u8; 8];
        buf[..rem.len()].copy_from_slice(rem);
        let w = u64::from_le_bytes(buf);
        h = (h.rotate_left(5) ^ w).wrapping_mul(K);
    }
    // final avalanche
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h
}

const EMPTY: u32 = u32::MAX;

pub struct KmerTable {
    k: usize,
    keys: Vec<u32>, // arena offset or EMPTY
    vals: Vec<u32>,
    mask: usize,
    len: usize,
    /// slots touched since the last clear (for O(len) clearing)
    touched: Vec<u32>,
}

impl KmerTable {
    pub fn new(k: usize, expected: usize) -> Self {
        let mut cap = 16usize;
        while cap < expected * 2 {
            cap <<= 1;
        }
        KmerTable { k, keys: vec![EMPTY; cap], vals: vec![0; cap], mask: cap - 1, len: 0, touched: Vec::new() }
    }

    pub fn clear(&mut self) {
        if self.touched.len() * 4 < self.keys.len() {
            for &i in &self.touched {
                self.keys[i as usize] = EMPTY;
            }
        } else {
            for x in self.keys.iter_mut() {
                *x = EMPTY;
            }
        }
        self.touched.clear();
        self.len = 0;
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    fn find_slot(&self, arena: &[u8], kmer: &[u8], h: u64) -> (usize, bool) {
        let mut i = (h as usize) & self.mask;
        loop {
            let key = self.keys[i];
            if key == EMPTY {
                return (i, false);
            }
            let off = key as usize;
            if &arena[off..off + self.k] == kmer {
                return (i, true);
            }
            i = (i + 1) & self.mask;
        }
    }

    pub fn get(&self, arena: &[u8], kmer: &[u8]) -> Option<u32> {
        let (i, found) = self.find_slot(arena, kmer, hash_kmer(kmer));
        if found {
            Some(self.vals[i])
        } else {
            None
        }
    }

    pub fn contains(&self, arena: &[u8], kmer: &[u8]) -> bool {
        self.find_slot(arena, kmer, hash_kmer(kmer)).1
    }

    /// Insert key (given by arena offset) if absent. Returns true if inserted.
    pub fn insert_if_absent(&mut self, arena: &[u8], off: u32, val: u32) -> bool {
        if (self.len + 1) * 2 > self.keys.len() {
            self.grow(arena);
        }
        let kmer = &arena[off as usize..off as usize + self.k];
        let (i, found) = self.find_slot(arena, kmer, hash_kmer(kmer));
        if found {
            return false;
        }
        self.keys[i] = off;
        self.vals[i] = val;
        self.len += 1;
        self.touched.push(i as u32);
        true
    }

    fn grow(&mut self, arena: &[u8]) {
        let new_cap = self.keys.len() * 2;
        let old_keys = std::mem::replace(&mut self.keys, vec![EMPTY; new_cap]);
        let old_vals = std::mem::replace(&mut self.vals, vec![0; new_cap]);
        self.mask = new_cap - 1;
        self.touched.clear();
        for (key, val) in old_keys.into_iter().zip(old_vals) {
            if key != EMPTY {
                let kmer = &arena[key as usize..key as usize + self.k];
                let (i, _) = self.find_slot(arena, kmer, hash_kmer(kmer));
                self.keys[i] = key;
                self.vals[i] = val;
                self.touched.push(i as u32);
            }
        }
    }

    pub fn remove(&mut self, arena: &[u8], kmer: &[u8]) -> bool {
        let (mut i, found) = self.find_slot(arena, kmer, hash_kmer(kmer));
        if !found {
            return false;
        }
        self.keys[i] = EMPTY;
        self.len -= 1;
        // backward shift deletion
        let mut j = i;
        loop {
            j = (j + 1) & self.mask;
            let key = self.keys[j];
            if key == EMPTY {
                break;
            }
            let kk = &arena[key as usize..key as usize + self.k];
            let home = (hash_kmer(kk) as usize) & self.mask;
            // can the entry at j move to i? yes if home is not cyclically in (i, j]
            let dist_home_j = j.wrapping_sub(home) & self.mask;
            let dist_i_j = j.wrapping_sub(i) & self.mask;
            if dist_home_j >= dist_i_j {
                self.keys[i] = key;
                self.vals[i] = self.vals[j];
                self.touched.push(i as u32);
                self.keys[j] = EMPTY;
                i = j;
            }
        }
        true
    }
}

/// 2-bit code of an ACGT base, or None for any other byte.
#[inline]
pub fn base_code(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

pub trait PackedKey: Copy + Eq + Default {
    fn hash64(self) -> u64;
    fn push(self, code: u8, mask: Self) -> Self;
    fn mask_for(k: usize) -> Self;
}

impl PackedKey for u64 {
    #[inline]
    fn hash64(self) -> u64 {
        let h = self.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        h ^ (h >> 29)
    }
    #[inline]
    fn push(self, code: u8, mask: Self) -> Self {
        ((self << 2) | code as u64) & mask
    }
    fn mask_for(k: usize) -> Self {
        if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        }
    }
}

impl PackedKey for u128 {
    #[inline]
    fn hash64(self) -> u64 {
        let lo = self as u64;
        let hi = (self >> 64) as u64;
        let h = (lo ^ hi.rotate_left(32).wrapping_mul(0xC2B2_AE3D_27D4_EB4F)).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        h ^ (h >> 29)
    }
    #[inline]
    fn push(self, code: u8, mask: Self) -> Self {
        ((self << 2) | code as u128) & mask
    }
    fn mask_for(k: usize) -> Self {
        if k >= 64 {
            u128::MAX
        } else {
            (1u128 << (2 * k)) - 1
        }
    }
}

/// Open-addressing set of packed k-mer codes with O(1) clear (generation stamps).
pub struct PackedSet<K: PackedKey> {
    keys: Vec<K>,
    gens: Vec<u32>,
    gen: u32,
    mask: usize,
    len: usize,
}

impl<K: PackedKey> PackedSet<K> {
    pub fn new(expected: usize) -> Self {
        let mut cap = 16usize;
        while cap < expected * 2 {
            cap <<= 1;
        }
        PackedSet { keys: vec![K::default(); cap], gens: vec![0; cap], gen: 1, mask: cap - 1, len: 0 }
    }
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn clear(&mut self) {
        self.gen = self.gen.wrapping_add(1);
        if self.gen == 0 {
            for g in self.gens.iter_mut() {
                *g = 0;
            }
            self.gen = 1;
        }
        self.len = 0;
    }
    #[inline]
    pub fn contains(&self, key: K) -> bool {
        let mut i = (key.hash64() as usize) & self.mask;
        loop {
            if self.gens[i] != self.gen {
                return false;
            }
            if self.keys[i] == key {
                return true;
            }
            i = (i + 1) & self.mask;
        }
    }
    /// returns true if newly inserted
    #[inline]
    pub fn insert(&mut self, key: K) -> bool {
        if (self.len + 1) * 2 > self.keys.len() {
            self.grow();
        }
        let mut i = (key.hash64() as usize) & self.mask;
        loop {
            if self.gens[i] != self.gen {
                self.gens[i] = self.gen;
                self.keys[i] = key;
                self.len += 1;
                return true;
            }
            if self.keys[i] == key {
                return false;
            }
            i = (i + 1) & self.mask;
        }
    }
    fn grow(&mut self) {
        let cap = self.keys.len() * 2;
        let old_keys = std::mem::replace(&mut self.keys, vec![K::default(); cap]);
        let old_gens = std::mem::replace(&mut self.gens, vec![0; cap]);
        let old_gen = self.gen;
        self.gen = 1;
        self.mask = cap - 1;
        self.len = 0;
        for (k, g) in old_keys.into_iter().zip(old_gens) {
            if g == old_gen {
                self.insert(k);
            }
        }
    }
}
