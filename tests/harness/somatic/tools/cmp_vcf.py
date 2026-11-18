#!/usr/bin/env python3
"""Compare two Mutect2 VCFs: header (ignoring the command-line line) and records field by field.
usage: cmp_vcf.py truth.vcf test.vcf [max_examples]"""
import sys
from collections import Counter

def load(p):
    hdr, recs = [], []
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('##'):
            if not (l.startswith('##GATKCommandLine') or l.startswith('##variant-rsCommandLine')):
                hdr.append(l)
        elif l.startswith('#'):
            hdr.append(l)
        else:
            recs.append(l.split('\t'))
    return hdr, recs

def info(s):
    d = {}
    for kv in s.split(';'):
        k, _, v = kv.partition('=')
        d[k] = v
    return d

a_h, a = load(sys.argv[1]); b_h, b = load(sys.argv[2])
mx = int(sys.argv[3]) if len(sys.argv) > 3 else 10
print(f"header: {'identical' if a_h == b_h else 'DIFFERS'}")
if a_h != b_h:
    for x in a_h:
        if x not in b_h: print('  only truth:', x[:150])
    for x in b_h:
        if x not in a_h: print('  only test :', x[:150])
ka = {(r[0], int(r[1]), r[3], r[4]): r for r in a}
kb = {(r[0], int(r[1]), r[3], r[4]): r for r in b}
only_a = sorted(set(ka) - set(kb)); only_b = sorted(set(kb) - set(ka))
diff = Counter(); ex = []
same = 0
for k in sorted(set(ka) & set(kb)):
    ra, rb = ka[k], kb[k]
    if ra == rb:
        same += 1; continue
    ia, ib = info(ra[7]), info(rb[7])
    for f in sorted(set(ia) | set(ib)):
        if ia.get(f) != ib.get(f):
            diff['INFO/' + f] += 1
            if len(ex) < mx: ex.append(f"{k[0]}:{k[1]} {f}: {ia.get(f)} vs {ib.get(f)}")
    if ra[8] != rb[8]:
        diff['FORMAT'] += 1
        if len(ex) < mx: ex.append(f"{k[0]}:{k[1]} FORMAT {ra[8]} vs {rb[8]}")
    else:
        keys = ra[8].split(':')
        for si in range(9, len(ra)):
            va, vb = ra[si].split(':'), rb[si].split(':')
            for j, f in enumerate(keys):
                x = va[j] if j < len(va) else None; y = vb[j] if j < len(vb) else None
                if x != y:
                    diff[f'S{si-9}/' + f] += 1
                    if len(ex) < mx: ex.append(f"{k[0]}:{k[1]} sample{si-9} {f}: {x} vs {y}")
    for j in (2, 5, 6):
        if ra[j] != rb[j]: diff[f'col{j}'] += 1
order_ok = [ (r[0], int(r[1]), r[3], r[4]) for r in a if (r[0], int(r[1]), r[3], r[4]) in kb] == [ (r[0], int(r[1]), r[3], r[4]) for r in b if (r[0], int(r[1]), r[3], r[4]) in ka]
print(f"records: truth {len(a)} test {len(b)}; identical lines {same}; only truth {len(only_a)}; only test {len(only_b)}; order {'same' if order_ok else 'DIFFERS'}")
for k in only_a[:mx]: print('  only truth:', '\t'.join(ka[k][:8])[:200])
for k in only_b[:mx]: print('  only test :', '\t'.join(kb[k][:8])[:200])
if diff: print('field differences:', dict(diff))
for e in ex: print('  ', e)
