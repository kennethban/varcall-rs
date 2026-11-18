#!/usr/bin/env python3
"""Synthetic test data for the self-contained tests: a small genome, a germline sample and a tumor/normal pair.

    simulate.py OUTDIR

Standard library only, seeded: the same files, byte for byte, on every machine (gzip with mtime 0). Writes
  ref.fa, ref.fa.fai, ref.dict          3 contigs, ~2 Mb: GC drift, STRs, homopolymers, an interspersed repeat
                                        family (multi-mapping reads, XA tags, MAPQ 0), segmental duplications
                                        (one across contigs), N gaps
  germline_L{1,2}_R{1,2}.fq.gz          sample GERM, two lanes, 2x150, ~25x
  tumor_R{1,2}.fq.gz, normal_R{1,2}.fq.gz   TUMOR (~40x; somatic subclones at AF 0.2 / 0.1 / 0.05) and NORMAL (~25x)
  known_sites.vcf                       BQSR known sites: most germline variants and some sites that are not
  germline_resource.vcf                 Mutect2 germline resource: germline variants with a population AF
  planted.tsv                           every planted variant (sample, contig, pos, ref, alt, genotype or AF)
Reads: Illumina-style names (tile, x, y), cycle-dependent qualities with low-quality tails, substitution errors
following the qualities and more after GG (context for BQSR), adapter read-through on short fragments, N calls,
PCR duplicates and optical duplicates (same tile, nearby x/y)."""
import bisect, gzip, hashlib, itertools, math, os, random, sys

SEED = 20261117
LINE = 60
BASES = 'ACGT'
COMP = str.maketrans('ACGTN', 'TGCAN')
ADAPTER = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG'  # TruSeq read 1 adapter


PERR = [10 ** (-q / 10) for q in range(64)]
QCHR = [chr(q + 33) for q in range(64)]
ALT = {b: [x for x in BASES if x != b] for b in BASES}
ALT['N'] = list(BASES)


def rc(s):
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------------------------------------------------------
# genome
# ---------------------------------------------------------------------------------------------------------------
def random_seq(rng, n, gc):
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return ''.join(rng.choices(BASES, weights=w, k=n))


def mutate(rng, s, rate):
    out = list(s)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in BASES if b != out[i]])
    return ''.join(out)


def make_genome(rng):
    sizes = [('chr1', 1_000_000), ('chr2', 600_000), ('chr3', 300_000)]
    contigs = {}
    for name, n in sizes:
        parts = []
        while sum(map(len, parts)) < n:
            parts.append(random_seq(rng, 1000, rng.uniform(0.35, 0.6)))
        contigs[name] = list(''.join(parts)[:n])
    # interspersed repeat family: ~150 copies of a 300 bp element, 2-15% diverged, either strand
    element = random_seq(rng, 300, 0.55)
    for _ in range(150):
        c = rng.choice(list(contigs))
        e = mutate(rng, element, rng.uniform(0.02, 0.15))
        if rng.random() < 0.5:
            e = rc(e)
        p = rng.randrange(5000, len(contigs[c]) - 5000)
        contigs[c][p:p + len(e)] = e
    # STRs (units 1-6 bp, 8-40 copies) every ~4 kb, and homopolymers
    for c, s in contigs.items():
        for p in range(3000, len(s) - 3000, 4000):
            p += rng.randrange(-1000, 1000)
            unit = random_seq(rng, rng.randint(1, 6), 0.5)
            rep = unit * rng.randint(8, 40)
            s[p:p + len(rep)] = rep
        for p in range(5000, len(s) - 5000, 25000):
            n = rng.randint(10, 15)
            s[p:p + n] = rng.choice(BASES) * n
    # segmental duplications: 10 kb, ~1% diverged; one within chr1, one chr1 -> chr2
    for (c1, p1, c2, p2) in [('chr1', 200_000, 'chr1', 700_000), ('chr1', 400_000, 'chr2', 300_000)]:
        seg = mutate(rng, ''.join(contigs[c1][p1:p1 + 10_000]), 0.01)
        contigs[c2][p2:p2 + 10_000] = seg
    # N gaps: chr1 middle (10 kb), chr2 start (2 kb, a telomere)
    contigs['chr1'][500_000:510_000] = 'N' * 10_000
    contigs['chr2'][0:2000] = 'N' * 2000
    return {c: ''.join(s) for c, s in contigs.items()}


def write_ref(out, contigs):
    with open(f'{out}/ref.fa', 'w') as f, open(f'{out}/ref.fa.fai', 'w') as fai, open(f'{out}/ref.dict', 'w') as d:
        d.write('@HD\tVN:1.6\n')
        off = 0
        for c, s in contigs.items():
            head = f'>{c}\n'
            f.write(head)
            off += len(head)
            fai.write(f'{c}\t{len(s)}\t{off}\t{LINE}\t{LINE + 1}\n')
            for i in range(0, len(s), LINE):
                f.write(s[i:i + LINE] + '\n')
            off += len(s) + (len(s) + LINE - 1) // LINE
            d.write(f'@SQ\tSN:{c}\tLN:{len(s)}\tM5:{hashlib.md5(s.encode()).hexdigest()}\n')


# ---------------------------------------------------------------------------------------------------------------
# variants and haplotypes
# ---------------------------------------------------------------------------------------------------------------
def str_loci(seq):
    """start, unit, end of tandem repeats of units 1-6 bp spanning >= 12 bp"""
    out = []
    i = 0
    while i < len(seq) - 12:
        best = None
        for u in range(1, 7):
            unit = seq[i:i + u]
            if 'N' in unit:
                break
            j = i + u
            while seq[j:j + u] == unit:
                j += u
            if j - i >= 12 and (best is None or j - i > best[2] - best[0]):
                best = (i, unit, j)
        if best:
            out.append(best)
            i = best[2]
        else:
            i += 1
    return out


def plant(rng, contigs, snp_rate, indel_rate, avoid):
    """Non-overlapping variants (contig, 0-based pos, ref, alt) on the reference, kept 20 bp apart; indels partly in
    STRs (unit expansions / contractions)"""
    vs = []
    for c, s in contigs.items():
        taken = set(avoid.get(c, ()))
        def free(p, n):
            return all(q not in taken for q in range(p - 20, p + n + 20)) and 'N' not in s[max(0, p - 1):p + n + 1]
        def take(p, n):
            taken.update(range(p - 20, p + n + 20))
        for p in range(1, len(s) - 1):
            if rng.random() < snp_rate and free(p, 1):
                r = s[p]
                if rng.random() < 0.05 and free(p, 3):  # MNP
                    ref = s[p:p + 2]
                    alt = ''.join(rng.choice([b for b in BASES if b != x]) for x in ref)
                    vs.append((c, p, ref, alt))
                    take(p, 2)
                else:
                    vs.append((c, p, r, rng.choice([b for b in BASES if b != r])))
                    take(p, 1)
            elif rng.random() < indel_rate and free(p, 25):
                n = min(int(rng.expovariate(0.35)) + 1, 20)
                if rng.random() < 0.5:   # deletion: ref = base before + deleted bases
                    vs.append((c, p, s[p:p + 1 + n], s[p]))
                else:
                    vs.append((c, p, s[p], s[p] + random_seq(rng, n, 0.5)))
                take(p, n + 1)
        # STR indels
        for (b, unit, e) in str_loci(s):
            if rng.random() < 0.25 and free(b, e - b):
                k = rng.randint(1, 3)
                if rng.random() < 0.5 and (e - b) // len(unit) > k:   # a contraction keeps at least one unit
                    vs.append((c, b - 1, s[b - 1:b + len(unit) * k], s[b - 1]))
                else:
                    vs.append((c, b - 1, s[b - 1], s[b - 1] + unit * k))
                take(b - 1, e - b + 1)
    vs.sort()
    return vs


def apply(seq, vs):
    out = []
    last = 0
    for (_, p, ref, alt) in vs:
        assert seq[p:p + len(ref)] == ref
        out.append(seq[last:p])
        out.append(alt)
        last = p + len(ref)
    out.append(seq[last:])
    return ''.join(out)


# ---------------------------------------------------------------------------------------------------------------
# reads
# ---------------------------------------------------------------------------------------------------------------
class Reads:
    """Illumina-like 2x150 pairs from a set of weighted haplotypes"""

    def __init__(self, rng, lane, flowcell):
        self.rng, self.lane, self.fc, self.n = rng, lane, flowcell, 0
        # cycle quality profiles: Q ~37 falling to ~30 over the read, binned like NovaSeq (2, 12, 23, 37) in some
        self.profiles = []
        for k in range(64):
            q = []
            for i in range(150):
                m = 37 - 7 * (i / 150) ** 2 - rng.random() * 3
                v = max(2, min(41, int(rng.gauss(m, 3))))
                q.append(v)
            if k % 4 == 0:
                q = [2 if v < 7 else 12 if v < 18 else 23 if v < 30 else 37 for v in q]
            self.profiles.append(q)
        self.cums = [list(itertools.accumulate(PERR[v] for v in q)) for q in self.profiles]

    def poisson(self, lam):
        lim, k, p = math.exp(-lam), 0, 1.0
        while True:
            p *= self.rng.random()
            if p <= lim:
                return k
            k += 1

    def read(self, frag):
        """(bases, qualities) of a 150-base read of `frag` (padded with adapter): substitution errors as the
        qualities say (Poisson count, positions by error weight), 3x as many right after GG, a low-quality tail in 4%
        of reads, N calls (quality 2) at 0.05%"""
        rng, n = self.rng, 150
        k = rng.randrange(len(self.profiles))
        q, cum = self.profiles[k], self.cums[k]
        if rng.random() < 0.04:
            t = rng.randint(5, 40)
            q = q[:n - t] + [2] * t
            cum = list(itertools.accumulate(PERR[v] for v in q))
        else:
            q = list(q)
        src = (frag + ADAPTER + 'A' * 150)[:n]
        s = list(src)
        lam = cum[-1]
        for _ in range(self.poisson(lam)):
            i = min(bisect.bisect(cum, rng.random() * lam), n - 1)
            s[i] = rng.choice(ALT[s[i]])
        j = src.find('GG')
        while j != -1 and j + 2 < n:
            i = j + 2
            if rng.random() < 2 * PERR[q[i]]:
                s[i] = rng.choice(ALT[s[i]])
            j = src.find('GG', j + 1)
        for _ in range(self.poisson(n * 0.0005)):
            i = rng.randrange(n)
            s[i] = 'N'
            q[i] = 2
        return ''.join(s), ''.join(QCHR[v] for v in q)

    def name(self, tile=None, xy=None):
        self.n += 1
        tile = tile or self.rng.choice([1101, 1102, 1201, 2101, 2202])
        x, y = xy or (self.rng.randint(1000, 30000), self.rng.randint(1000, 30000))
        return f'SIM:1:{self.fc}:{self.lane}:{tile}:{x}:{y}', tile, (x, y)

    def pairs(self, haps, depth, out1, out2):
        """haps: [(weight, {contig: seq})]; depth over the summed contig length"""
        rng = self.rng
        total = sum(len(s) for s in haps[0][1].values())
        n_frag = int(depth * total / 300)
        weights = [w for w, _ in haps]
        contig_names = list(haps[0][1])
        lens = [len(haps[0][1][c]) for c in contig_names]
        with gzip.GzipFile(out1, 'wb', mtime=0) as f1, gzip.GzipFile(out2, 'wb', mtime=0) as f2:
            for _ in range(n_frag):
                hap = rng.choices(haps, weights=weights)[0][1]
                c = rng.choices(contig_names, weights=lens)[0]
                seq = hap[c]
                flen = max(80, min(800, int(rng.gauss(350, 70))))
                p = rng.randrange(0, max(1, len(seq) - flen))
                frag = seq[p:p + flen]
                if frag.count('N') > 20:
                    continue
                if rng.random() < 0.5:
                    frag = rc(frag)
                copies = 1
                r = rng.random()
                if r < 0.04:
                    copies = 2          # PCR duplicate
                elif r < 0.05:
                    copies = 3
                optical = rng.random() < 0.01
                tile = xy = None
                for k in range(copies + (1 if optical else 0)):
                    if k > 0 and optical and k == copies:
                        xy = (xy[0] + rng.randint(-40, 40), xy[1] + rng.randint(-40, 40))
                        nm, _, _ = self.name(tile, xy)
                    else:
                        nm, tile, xy = self.name()
                    s1, q1 = self.read(frag)
                    s2, q2 = self.read(rc(frag))
                    for f, s, q, m in ((f1, s1, q1, 1), (f2, s2, q2, 2)):
                        f.write(f'@{nm} {m}:N:0:ACGTACGT\n{s}\n+\n{q}\n'.encode())


# ---------------------------------------------------------------------------------------------------------------
def write_vcf(path, contigs, rows, info_header=()):
    with open(path, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        for h in info_header:
            f.write(h + '\n')
        for c, s in contigs.items():
            f.write(f'##contig=<ID={c},length={len(s)}>\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        order = {c: i for i, c in enumerate(contigs)}
        for (c, p, ref, alt, info) in sorted(rows, key=lambda r: (order[r[0]], r[1])):
            f.write(f'{c}\t{p + 1}\t.\t{ref}\t{alt}\t.\t.\t{info}\n')


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    rng = random.Random(SEED)
    contigs = make_genome(rng)
    write_ref(out, contigs)

    germ = plant(rng, contigs, 1 / 1000, 1 / 8000, {})
    # genotypes: 40% hom alt, else het on haplotype A or B; ~3% of het sites get a second alt on the other haplotype
    hap_a, hap_b, planted = [], [], []
    multi_b = []
    for v in germ:
        r = rng.random()
        if r < 0.4:
            hap_a.append(v); hap_b.append(v); planted.append(('germline',) + v + ('1/1',))
        elif r < 0.7:
            hap_a.append(v); planted.append(('germline',) + v + ('0/1',))
            if len(v[2]) == 1 and len(v[3]) == 1 and rng.random() < 0.1:
                alt2 = rng.choice([b for b in BASES if b not in (v[2], v[3])])
                w = (v[0], v[1], v[2], alt2)
                multi_b.append(w); planted.append(('germline',) + w + ('1/2',))
        else:
            hap_b.append(v); planted.append(('germline',) + v + ('0/1',))
    hap_b = sorted(hap_b + multi_b)
    germ_avoid = {}
    for (c, p, ref, alt) in germ:
        germ_avoid.setdefault(c, []).extend(range(p, p + len(ref)))
    soma = plant(rng, contigs, 1 / 20000, 1 / 150000, germ_avoid)
    clones = [[], [], []]
    for v in soma:
        k = rng.randrange(3)
        clones[k].append(v)
        planted.append(('somatic',) + v + (('0.2', '0.1', '0.05')[k],))

    def hap(vs):
        by = {c: [] for c in contigs}
        for v in vs:
            by[v[0]].append(v)
        return {c: apply(contigs[c], sorted(by[c])) for c in contigs}

    A, B = hap(hap_a), hap(hap_b)
    A1, A2, A3 = hap(sorted(hap_a + clones[0])), hap(sorted(hap_a + clones[1])), hap(sorted(hap_a + clones[2]))

    Reads(rng, 1, 'FCGERM1').pairs([(1, A), (1, B)], 12.5, f'{out}/germline_L1_R1.fq.gz', f'{out}/germline_L1_R2.fq.gz')
    Reads(rng, 2, 'FCGERM1').pairs([(1, A), (1, B)], 12.5, f'{out}/germline_L2_R1.fq.gz', f'{out}/germline_L2_R2.fq.gz')
    Reads(rng, 1, 'FCNORM1').pairs([(1, A), (1, B)], 25, f'{out}/normal_R1.fq.gz', f'{out}/normal_R2.fq.gz')
    # tumor: reads from hap B half the time; hap A splits into the clones (somatic AF 0.2, 0.1, 0.05) and plain A
    Reads(rng, 1, 'FCTUMR1').pairs([(0.5, B), (0.2, A1), (0.1, A2), (0.05, A3), (0.15, A)], 40,
                                    f'{out}/tumor_R1.fq.gz', f'{out}/tumor_R2.fq.gz')

    known = [(c, p, r, a, '.') for (c, p, r, a) in germ if rng.random() < 0.7]
    for _ in range(300):  # known sites that are not variant in these samples
        c = rng.choice(list(contigs))
        p = rng.randrange(3000, len(contigs[c]) - 3000)
        b = contigs[c][p]
        if b != 'N' and all(abs(p - q) > 50 for (cc, q, _, _) in germ if cc == c):
            known.append((c, p, b, rng.choice([x for x in BASES if x != b]), '.'))
    dedup = {}
    for row in known:
        dedup.setdefault((row[0], row[1]), row)
    write_vcf(f'{out}/known_sites.vcf', contigs, list(dedup.values()))
    res = [(c, p, r, a, f'AF={rng.choice([0.0005, 0.002, 0.01, 0.05, 0.2, 0.5])}') for (c, p, r, a) in germ]
    write_vcf(f'{out}/germline_resource.vcf', contigs, res,
              ['##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'])
    with open(f'{out}/planted.tsv', 'w') as f:
        f.write('sample\tcontig\tpos\tref\talt\tgenotype_or_af\n')
        for (s, c, p, r, a, g) in planted:
            f.write(f'{s}\t{c}\t{p + 1}\t{r}\t{a}\t{g}\n')


if __name__ == '__main__':
    main()
