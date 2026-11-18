//! GATK's global `java.util.Random` (`Utils.getRandomGenerator()`, seed 47382911) for one *domain* (one
//! GATK JVM: the whole run, or one shard file). Its two users, positional downsampling (`nextInt` per read
//! beyond 50 at a start position) and the QD jitter (`nextGaussian`, VCF mode), draw from it in GATK's
//! traversal order: contigs in dictionary order, and inside a contig in the order fixed by
//! `hc_downsample`'s schedule. The generator is a baton handed from contig to contig of the domain.

use hc_downsample::JavaRandom;
use std::sync::{Condvar, Mutex};

struct State<T> {
    /// index (within the domain) of the contig that may use the generator now
    cursor: usize,
    /// None while a contig has it checked out
    rng: Option<T>,
    /// contigs that never use the generator (no drawing group, no QD jitter)
    pass: Vec<bool>,
}

/// A generator handed from contig to contig of a domain. `T`: GATK's `java.util.Random` (default) or, for
/// contamination downsampling, its commons-math `Well19937c` (`Utils.getRandomDataGenerator()`), which no other
/// HaplotypeCaller step draws from.
pub struct DomainRng<T = JavaRandom> {
    m: Mutex<State<T>>,
    cv: Condvar,
}

impl DomainRng<JavaRandom> {
    pub fn new(n_contigs: usize) -> DomainRng {
        DomainRng::with(n_contigs, JavaRandom::gatk())
    }
}

impl<T> DomainRng<T> {
    pub fn with(n_contigs: usize, init: T) -> DomainRng<T> {
        DomainRng { m: Mutex::new(State { cursor: 0, rng: Some(init), pass: vec![false; n_contigs] }), cv: Condvar::new() }
    }

    fn skip(s: &mut State<T>) {
        while s.cursor < s.pass.len() && s.pass[s.cursor] && s.rng.is_some() {
            s.cursor += 1;
        }
    }

    /// Blocks until contig `k`'s turn (every earlier contig of the domain is done with the generator).
    pub fn take(&self, k: usize) -> T {
        let mut s = self.m.lock().unwrap();
        loop {
            if s.cursor == k && s.rng.is_some() {
                return s.rng.take().unwrap();
            }
            s = self.cv.wait(s).unwrap();
        }
    }

    pub fn try_take(&self, k: usize) -> Option<T> {
        let mut s = self.m.lock().unwrap();
        if s.cursor == k && s.rng.is_some() {
            s.rng.take()
        } else {
            None
        }
    }

    /// Contig `k` is done with the generator.
    pub fn finish(&self, k: usize, rng: T) {
        let mut s = self.m.lock().unwrap();
        assert!(s.cursor == k && s.rng.is_none(), "RNG baton out of order");
        s.rng = Some(rng);
        s.cursor = k + 1;
        Self::skip(&mut s);
        self.cv.notify_all();
    }

    /// Contig `k` makes no draw at all.
    pub fn pass(&self, k: usize) {
        let mut s = self.m.lock().unwrap();
        s.pass[k] = true;
        Self::skip(&mut s);
        self.cv.notify_all();
    }
}

/// Replaces every QD jitter placeholder (`hc_genotype::QD_JITTER_PLACEHOLDER`, byte 1) in `text`, in order,
/// with `QualByDepth.fixTooHighQD`'s value `30 + 3 * nextGaussian()` formatted `%.2f`. Returns the number of draws.
pub fn resolve_qd(text: &[u8], rng: &mut JavaRandom, out: &mut Vec<u8>) -> u64 {
    let mut n = 0;
    let mut last = 0;
    for i in memchr::memchr_iter(1u8, text) {
        out.extend_from_slice(&text[last..i]);
        out.extend_from_slice(hc_genotype::javamath::format_fixed(rng.qd_jitter(), 2).as_bytes());
        last = i + 1;
        n += 1;
    }
    out.extend_from_slice(&text[last..]);
    n
}
