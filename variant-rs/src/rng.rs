//! GATK's global `java.util.Random` (`Utils.getRandomGenerator()`, seed 47382911) for one *domain* (one
//! GATK JVM: the whole run, or one shard file). Its two users, positional downsampling (`nextInt` per read
//! beyond 50 at a start position) and the QD jitter (`nextGaussian`, VCF mode), draw from it in GATK's
//! traversal order: contigs in dictionary order, and inside a contig in the order fixed by
//! `hc_downsample`'s schedule. The generator is a baton handed from contig to contig of the domain (`DomainRng`);
//! in GVCF mode each contig instead gets its starting state from a replay of the earlier contigs' draws (`Replay`).

use hc_downsample::JavaRandom;
use std::sync::{Condvar, Mutex};

struct State<T> {
    /// index (within the domain) of the contig that may use the generator now
    cursor: usize,
    /// None while a contig has it checked out
    rng: Option<T>,
    /// contigs that never use the generator (no drawing group, no QD jitter)
    pass: Vec<bool>,
    /// segments whose range an earlier segment took over (see `set_absorbs`): they never get the generator
    absorbed: Vec<bool>,
    /// number of following segments each segment took over, and how many of them are marked `absorbed`
    absorbs: Vec<usize>,
    applied: Vec<usize>,
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
        DomainRng { m: Mutex::new(State { cursor: 0, rng: Some(init), pass: vec![false; n_contigs], absorbed: vec![false; n_contigs], absorbs: vec![0; n_contigs], applied: vec![0; n_contigs] }), cv: Condvar::new() }
    }

    /// The segments `k` took over become void once it is `k`'s turn (then `k` itself is not void: a void
    /// segment's own decisions never apply).
    fn apply(s: &mut State<T>, k: usize) {
        while s.applied[k] < s.absorbs[k] {
            s.applied[k] += 1;
            let j = k + s.applied[k];
            s.absorbed[j] = true;
            s.pass[j] = true;
        }
    }

    fn skip(s: &mut State<T>) {
        while s.cursor < s.pass.len() {
            let c = s.cursor;
            if !s.absorbed[c] {
                Self::apply(s, c);
            }
            if s.pass[c] && s.rng.is_some() {
                s.cursor += 1;
            } else {
                break;
            }
        }
    }

    /// Blocks until contig `k`'s turn (every earlier contig of the domain is done with the generator).
    /// `None`: `k` is a segment the previous one took over (`absorb`); its results are void.
    pub fn take(&self, k: usize) -> Option<T> {
        let mut s = self.m.lock().unwrap();
        loop {
            if s.absorbed[k] {
                return None;
            }
            if s.cursor == k && s.rng.is_some() {
                return s.rng.take();
            }
            s = self.cv.wait(s).unwrap();
        }
    }

    pub fn try_take(&self, k: usize) -> Option<T> {
        let mut s = self.m.lock().unwrap();
        if !s.absorbed[k] && s.cursor == k && s.rng.is_some() {
            s.rng.take()
        } else {
            None
        }
    }

    /// Contig `k` is done with the generator.
    pub fn finish(&self, k: usize, rng: T) {
        let mut s = self.m.lock().unwrap();
        assert!(s.cursor == k && s.rng.is_none(), "RNG baton out of order");
        Self::apply(&mut s, k);
        s.rng = Some(rng);
        s.cursor = k + 1;
        Self::skip(&mut s);
        self.cv.notify_all();
    }

    /// Segment `k` processes the ranges of the next `n` segments: they are skipped and never get the generator
    /// (their results are void). Applied when it is `k`'s turn, so that a segment that is itself void decides nothing.
    pub fn set_absorbs(&self, k: usize, n: usize) {
        let mut s = self.m.lock().unwrap();
        s.absorbs[k] = n;
        if s.cursor == k && !s.absorbed[k] {
            Self::apply(&mut s, k);
        }
        Self::skip(&mut s);
        self.cv.notify_all();
    }

    pub fn is_absorbed(&self, k: usize) -> bool {
        self.m.lock().unwrap().absorbed[k]
    }

    /// Index of the contig whose turn it is (every earlier one is done with the generator).
    pub fn turn(&self) -> usize {
        self.m.lock().unwrap().cursor
    }

    /// Waits up to `ms` milliseconds for the turn to move.
    pub fn wait_turn_change(&self, seen: usize, ms: u64) {
        let s = self.m.lock().unwrap();
        if s.cursor == seen {
            let _ = self.cv.wait_timeout(s, std::time::Duration::from_millis(ms)).unwrap();
        }
    }

    /// Contig `k` makes no draw at all.
    pub fn pass(&self, k: usize) {
        let mut s = self.m.lock().unwrap();
        s.pass[k] = true;
        Self::skip(&mut s);
        self.cv.notify_all();
    }
}

/// GVCF mode: the generator's only users are the downsampler's drawing groups (no QD jitter, no contamination), and
/// how many reads a group has does not depend on the generator. So the state each contig of a domain starts with
/// follows from the group sizes of the contigs before it: a pre-scan (`frontend::count_drawing_groups`, in chunks,
/// in parallel) collects them, the replay makes their draws in order (`reservoir_draws`), and every contig gets its
/// starting state as soon as the contigs before it are scanned, instead of when they are done. At its end each
/// contig compares its generator with the replay's state for the next contig (`check`); a difference stops the run.
/// The last contig of a domain is not scanned (no contig starts after it).
pub struct Replay {
    target: usize,
    t0: std::time::Instant,
    m: Mutex<ReplayState>,
    cv: Condvar,
}

struct ReplayState {
    /// per contig (all but the last), per chunk: the chunk's drawing groups (start, reads) once scanned
    chunks: Vec<Vec<Option<Vec<(i64, u32)>>>>,
    /// per contig: the generator's state when it starts
    starts: Vec<Option<JavaRandom>>,
    rng: JavaRandom,
    /// next contig to replay
    next: usize,
    groups: u64,
    draws: u64,
    failed: bool,
}

impl Replay {
    /// `chunks[k]`: number of scan chunks of contig k, for the domain's contigs but the last.
    pub fn new(chunks: &[usize], target: usize, t0: std::time::Instant) -> Replay {
        let mut starts = vec![None; chunks.len() + 1];
        starts[0] = Some(JavaRandom::gatk());
        let st = ReplayState { chunks: chunks.iter().map(|&n| vec![None; n]).collect(), starts, rng: JavaRandom::gatk(), next: 0, groups: 0, draws: 0, failed: false };
        let r = Replay { target, t0, m: Mutex::new(st), cv: Condvar::new() };
        r.advance(&mut r.m.lock().unwrap());
        r
    }

    fn advance(&self, s: &mut ReplayState) {
        while s.next < s.chunks.len() && s.chunks[s.next].iter().all(|c| c.is_some()) {
            let k = s.next;
            for c in std::mem::take(&mut s.chunks[k]) {
                for (_, n) in c.unwrap() {
                    s.groups += 1;
                    s.draws += n as u64 - self.target as u64;
                    hc_downsample::reservoir_draws(n as usize, self.target, &mut s.rng);
                }
            }
            s.next = k + 1;
            s.starts[k + 1] = Some(s.rng.clone());
            if s.next == s.chunks.len() && s.next > 0 {
                eprintln!("replay: {} contig(s), {} drawing groups, {} draws, ready at {:.1}s", s.next + 1, s.groups, s.draws, self.t0.elapsed().as_secs_f64());
            }
        }
    }

    /// Chunk `c` of contig `k` is scanned: its drawing groups.
    pub fn add(&self, k: usize, c: usize, groups: Vec<(i64, u32)>) {
        let mut s = self.m.lock().unwrap();
        s.chunks[k][c] = Some(groups);
        self.advance(&mut s);
        self.cv.notify_all();
    }

    /// A scan failed: contigs waiting for a state stop.
    pub fn fail(&self) {
        self.m.lock().unwrap().failed = true;
        self.cv.notify_all();
    }

    fn wait(&self, k: usize) -> JavaRandom {
        let mut s = self.m.lock().unwrap();
        loop {
            if let Some(r) = &s.starts[k] {
                return r.clone();
            }
            assert!(!s.failed, "GVCF generator replay: the pre-scan failed");
            s = self.cv.wait(s).unwrap();
        }
    }

    /// The generator contig `k` starts with (blocks until the contigs before it are scanned).
    pub fn start(&self, k: usize) -> JavaRandom {
        self.wait(k)
    }

    /// Contig `k` is done: `end` is its generator (None: it never drew). Err if it differs from the replay's state
    /// for contig `k + 1`.
    pub fn check(&self, k: usize, end: Option<&JavaRandom>) -> Result<(), String> {
        if k + 1 >= self.m.lock().unwrap().starts.len() {
            return Ok(());
        }
        let want = self.wait(k + 1).state();
        let got = match end {
            Some(r) => r.state(),
            None => self.wait(k).state(),
        };
        if got == want {
            Ok(())
        } else {
            Err(format!("GVCF generator replay: contig {k} of the domain ends with the generator at {got:?}, the replay has {want:?} (VARIANT_RS_REPLAY=0 turns the replay off)"))
        }
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

#[cfg(test)]
mod tests {
    use super::{DomainRng, Replay};
    use hc_downsample::JavaRandom;

    #[test]
    fn replay_gives_each_contig_the_state_after_the_earlier_ones() {
        let r = Replay::new(&[2, 1], 50, std::time::Instant::now());
        let mut g = JavaRandom::gatk();
        assert_eq!(r.start(0).state(), g.state());
        // contig 1's state needs both chunks of contig 0
        r.add(0, 1, vec![(900, 70)]);
        r.add(1, 0, vec![(5, 51)]);
        r.add(0, 0, vec![(10, 60), (20, 50 + 1)]);
        for n in [60, 51, 70] {
            hc_downsample::reservoir_draws(n, 50, &mut g);
        }
        assert_eq!(r.start(1).state(), g.state());
        assert!(r.check(0, Some(&g)).is_ok());
        assert!(r.check(0, None).is_err());
        let mut g2 = g.clone();
        hc_downsample::reservoir_draws(51, 50, &mut g2);
        assert_eq!(r.start(2).state(), g2.state());
        assert!(r.check(1, Some(&g)).is_err());
        assert!(r.check(2, None).is_ok()); // the last contig is not checked
    }

    #[test]
    fn absorbed_segments_are_skipped() {
        let d: DomainRng<u32> = DomainRng::with(4, 7);
        // segment 1 takes over segment 2 before its turn: applied on its turn
        d.set_absorbs(1, 1);
        assert!(!d.is_absorbed(2));
        let g = d.take(0).unwrap();
        d.finish(0, g + 1);
        assert!(d.is_absorbed(2));
        let g = d.take(1).unwrap();
        assert_eq!(g, 8);
        assert!(d.take(2).is_none());
        d.finish(1, g + 1);
        assert_eq!(d.turn(), 3);
        assert_eq!(d.take(3), Some(9));
    }

    #[test]
    fn void_segment_decides_nothing() {
        let d: DomainRng<u32> = DomainRng::with(4, 0);
        d.set_absorbs(0, 1); // 0 takes 1 over
        d.set_absorbs(1, 1); // 1 (void) would take 2 over: ignored
        assert!(d.is_absorbed(1));
        let g = d.take(0).unwrap();
        d.finish(0, g);
        assert!(!d.is_absorbed(2));
        assert_eq!(d.turn(), 2);
        d.pass(2);
        assert_eq!(d.turn(), 3);
    }
}
