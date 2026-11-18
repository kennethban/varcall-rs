#!/usr/bin/env python3
"""blocks.py A.bam B.bam: BGZF blocks of two files compared as raw bytes: count of blocks and of differing ones."""
import sys, struct
def blocks(p):
    d = open(p, 'rb').read(); o = 0; out = []
    while o < len(d):
        bs = struct.unpack('<H', d[o+16:o+18])[0] + 1; out.append(d[o:o+bs]); o += bs
    return out
a, b = blocks(sys.argv[1]), blocks(sys.argv[2])
diff = [i for i in range(min(len(a), len(b))) if a[i] != b[i]]
print(f"{len(a)}/{len(b)} blocks, {len(diff)} differ" + (f" (blocks {diff[:8]})" if diff else ""))
