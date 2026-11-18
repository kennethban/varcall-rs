#!/usr/bin/env python3
"""Synthetic read sets for edge cases the real-read tests do not reach (adversarial review, REVIEW_FINDINGS.md).
   synth.py SET OUT_PREFIX [N]   -> OUT_PREFIX.R1.fq.gz, OUT_PREFIX.R2.fq.gz (bwa-mem2 gives the truth)
Sets:
  longdel   2x250 pairs; R1 = 150 bp of reference, a 75-100 bp deletion, 100 bp more (the right extension's second
            band, finding 2), either strand; R2 a plain mate 400-650 bp downstream.
  unmapped  random 2x250 pairs with 200-character names: nothing maps, so every record lands in the BAM writer's
            unplaced bucket (finding 3: > 4 GiB in one bucket). N pairs default 3.6M (~4.3 GiB of BAM records).
  smallref  2x150 pairs from tests/data/small.fa (chr20:10-15 Mb as its own reference, bwa-mem2 index; a reference
            under ~500 Mb crashed --gpu start-up), with empty mates (zero-length reads crashed mate rescue).
  longread  pairs of 700-3000 bp reads (re-scoring of short seeds from 725 bp on: crashed before 2026-11-17).
Reference: ref_hg38 (Broad hg38), positions from a fixed seed on chr1-chr22, windows without N."""
import gzip, os, random, subprocess, sys

REF = os.path.join(os.environ['VARCALL_DATA'], 'ref_hg38/Homo_sapiens_assembly38.fasta')
COMP = str.maketrans('ACGTN', 'TGCAN')

def fetch(chrom, beg, end):  # 0-based half-open
    out = subprocess.run(['samtools', 'faidx', REF, f'{chrom}:{beg + 1}-{end}'], capture_output=True, text=True, check=True).stdout
    return ''.join(out.split('\n')[1:]).upper()

def rc(s):
    return s.translate(COMP)[::-1]

def lengths():
    fai = {}
    for l in open(REF + '.fai'):
        c = l.split('\t')
        fai[c[0]] = int(c[1])
    return fai

def longdel(rng, fai, n):
    pairs = []
    chroms = [f'chr{i}' for i in range(1, 23)]
    while len(pairs) < n:
        c = rng.choice(chroms)
        p = rng.randrange(1_000_000, fai[c] - 2_000_000)
        d = rng.randrange(75, 101)
        win = fetch(c, p, p + 1000)
        if 'N' in win:
            continue
        r1 = win[0:150] + win[150 + d:250 + d]
        r2 = rc(win[400 + d:650 + d])
        name = f'{c}_{p}_del{d}'
        if rng.random() < 0.5:  # the pair from the other strand
            r1, r2 = rc(win[400 + d:650 + d]), rc(win[0:150] + win[150 + d:250 + d])
            name += '_rev'
        pairs.append((name, r1, r2))
    return pairs

def unmapped(rng, fai, n):
    import os
    lut = bytes(b'ACGT'[i & 3] for i in range(256))
    for i in range(n):
        b = os.urandom(500).translate(lut).decode()
        yield (f'u{i:09d}_' + 'x' * 190, b[:250], b[250:])

def smallref_reads(rng, fai, n, ref=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'small.fa')):
    """2x150 pairs from the small reference (tests/data/small.fa), 1% substitutions; every 50th pair has an
    empty mate (zero-length sequence: mate rescue with an empty query), every 47th an empty read 1"""
    seq = ''.join(l.strip() for l in open(ref) if not l.startswith('>')).upper()
    out = []
    for i in range(n):
        p = rng.randrange(0, len(seq) - 700)
        frag = seq[p:p + rng.randrange(300, 600)]
        if 'N' in frag:
            continue
        mut = lambda s: ''.join(rng.choice('ACGT') if rng.random() < 0.01 else c for c in s)
        r1, r2 = mut(frag[:150]), mut(rc(frag[-150:]))
        if rng.random() < 0.5:
            r1, r2 = r2, r1
        if i % 50 == 0:
            r2 = ''
        if i % 47 == 0:
            r1 = ''
        out.append((f'small_{p}', r1, r2))
    return out

def longread(rng, fai, n):
    """pairs of 700-3000 bp reads (mem_seed_sw from 725 bp on), 1% substitutions, a 1-5 bp indel per ~300 bp;
    every 10th pair is 724 or 725 bp, either side of the cutoff"""
    chroms = [f'chr{i}' for i in range(1, 23)]
    def mutate(s):
        out, i = [], 0
        while i < len(s):
            r = rng.random()
            if r < 0.01:
                out.append(rng.choice('ACGT')); i += 1
            elif r < 0.01 + 1 / 600:
                i += rng.randrange(1, 6)  # deletion
            elif r < 0.01 + 1 / 300:
                out.append(''.join(rng.choice('ACGT') for _ in range(rng.randrange(1, 6)))); out.append(s[i]); i += 1
            else:
                out.append(s[i]); i += 1
        return ''.join(out)
    pairs = []
    while len(pairs) < n:
        c = rng.choice(chroms)
        rl = rng.choice((724, 725)) if len(pairs) % 10 == 0 else rng.randrange(700, 3001)
        fl = rl * 2 + rng.randrange(0, 2000)
        p = rng.randrange(1_000_000, fai[c] - 2_000_000)
        win = fetch(c, p, p + fl)
        if 'N' in win:
            continue
        r1, r2 = mutate(win[:rl])[:rl], mutate(rc(win[-rl:]))[:rl]
        if rng.random() < 0.5:
            r1, r2 = r2, r1
        pairs.append((f'{c}_{p}_len{rl}', r1, r2))
    return pairs

SETS = {'longdel': longdel, 'unmapped': unmapped, 'smallref': smallref_reads, 'longread': longread}

def main():
    which, prefix = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    pairs = SETS[which](random.Random(which), lengths(), n)
    f1 = gzip.open(f'{prefix}.R1.fq.gz', 'wt', compresslevel=1)
    f2 = gzip.open(f'{prefix}.R2.fq.gz', 'wt', compresslevel=1)
    for i, p in enumerate(pairs):
        for k, f in ((1, f1), (2, f2)):
            s = p[k]
            f.write(f'@syn{i}_{p[0]}/{k}\n{s}\n+\n{"I" * len(s)}\n')
    f1.close()
    f2.close()

main()
