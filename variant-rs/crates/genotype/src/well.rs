//! GATK's second global generator, `Utils.getRandomDataGenerator()`: commons-math 3.5
//! `RandomDataGenerator(new Well19937c(47382911))`. HaplotypeCaller draws from it only in contamination
//! downsampling (`AlleleBiasedDownsamplingUtils.downsampleElements` -> `MathUtils.sampleIndicesWithoutReplacement`
//! -> `RandomDataGenerator.nextPermutation`).

const K: usize = 19937;
const M1: usize = 70;
const M2: usize = 179;
const M3: usize = 449;
const R: usize = (K + 31) / 32;

#[derive(Clone)]
pub struct Well19937c {
    v: Vec<i32>,
    index: usize,
}

impl Well19937c {
    /// `new Well19937c(long seed)` (`AbstractWell.setSeed(int[])` from the two halves of the seed).
    pub fn new(seed: i64) -> Well19937c {
        let mut v = vec![0i32; R];
        let s = [(seed as u64 >> 32) as i32, (seed as u64 & 0xffff_ffff) as i32];
        v[..2].copy_from_slice(&s);
        for i in 2..R {
            let l = v[i - 2] as i64; // sign-extended, as Java's long l = v[i - seed.length]
            v[i] = ((1812433253i64.wrapping_mul(l ^ (l >> 30)).wrapping_add(i as i64)) & 0xffff_ffff) as i32;
        }
        Well19937c { v, index: 0 }
    }

    /// GATK's seed.
    pub fn gatk() -> Well19937c {
        Well19937c::new(47382911)
    }

    fn next(&mut self, bits: u32) -> i32 {
        let idx = self.index;
        let rm1 = (idx + R - 1) % R;
        let rm2 = (idx + R - 2) % R;
        let v = &mut self.v;
        let v0 = v[idx] as u32;
        let vm1 = v[(idx + M1) % R] as u32;
        let vm2 = v[(idx + M2) % R] as u32;
        let vm3 = v[(idx + M3) % R] as u32;
        let z0 = (0x8000_0000 & v[rm1] as u32) ^ (0x7FFF_FFFF & v[rm2] as u32);
        let z1 = (v0 ^ (v0 << 25)) ^ (vm1 ^ (vm1 >> 27));
        let z2 = (vm2 >> 9) ^ (vm3 ^ (vm3 >> 1));
        let z3 = z1 ^ z2;
        let mut z4 = z0 ^ (z1 ^ (z1 << 9)) ^ (z2 ^ (z2 << 21)) ^ (z3 ^ (z3 >> 21));
        v[idx] = z3 as i32;
        v[rm1] = z4 as i32;
        v[rm2] = (v[rm2] as u32 & 0x8000_0000) as i32;
        self.index = rm1;
        z4 ^= (z4 << 7) & 0xe46e_1700;
        z4 ^= (z4 << 15) & 0x9b86_8000;
        (z4 >> (32 - bits)) as i32
    }

    /// `BitsStreamGenerator.nextInt(int n)`
    pub fn next_int(&mut self, n: i32) -> i32 {
        assert!(n > 0);
        if (n & -n) == n {
            return ((n as i64 * self.next(31) as i64) >> 31) as i32;
        }
        loop {
            let bits = self.next(31);
            let val = bits % n;
            if bits.wrapping_sub(val).wrapping_add(n - 1) >= 0 {
                return val;
            }
        }
    }

    /// `RandomDataGenerator.nextPermutation(n, k)`: `MathArrays.shuffle(natural(n), rng)` (TAIL), first k.
    pub fn next_permutation(&mut self, n: usize, k: usize) -> Vec<usize> {
        assert!(k <= n && k > 0);
        let mut list: Vec<usize> = (0..n).collect();
        for i in (0..n).rev() {
            // UniformIntegerDistribution(rng, 0, i).sample() = 0 + nextInt(i + 1); i == 0 swaps with itself
            let target = if i == 0 { 0 } else { self.next_int(i as i32 + 1) as usize };
            list.swap(target, i);
        }
        list.truncate(k);
        list
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Values from commons-math 3.5 (jshell: `new Well19937c(47382911L)`, `nextInt(100)` x5, then
    /// `new RandomDataGenerator(new Well19937c(47382911L)).nextPermutation(10, 4)`); see well_check.jsh.
    #[test]
    fn matches_commons_math() {
        let expect: Vec<i32> = include_str!("../tests/well19937c.txt").split_whitespace().map(|x| x.parse().unwrap()).collect();
        let mut w = Well19937c::gatk();
        let ints: Vec<i32> = (0..5).map(|_| w.next_int(100)).collect();
        assert_eq!(ints, expect[..5]);
        let mut w = Well19937c::gatk();
        let p: Vec<i32> = w.next_permutation(10, 4).into_iter().map(|x| x as i32).collect();
        assert_eq!(p, expect[5..9]);
        let mut w = Well19937c::gatk();
        let big: Vec<i32> = (0..3).map(|_| w.next_int(1_000_000_007)).collect();
        assert_eq!(big, expect[9..12]);
    }
}
