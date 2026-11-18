#!/usr/bin/env python3
"""cmpbam.py A.bam B.bam: compare decompressed BAM records byte for byte (headers excluded)."""
import sys, struct, zlib
def records(path):
    f = open(path, 'rb'); buf = b''; pos = 0
    def blocks():
        while True:
            h = f.read(18)
            if len(h) < 18: return
            bsize = struct.unpack('<H', h[16:18])[0] + 1
            rest = f.read(bsize - 18)
            yield zlib.decompress(rest[:-8], -15)
    data = bytearray(); it = blocks()
    def need(n):
        while len(data) < n:
            try: data.extend(next(it))
            except StopIteration: return False
        return True
    need(8); lt = struct.unpack('<i', data[4:8])[0]; need(12 + lt)
    nref = struct.unpack('<i', data[8+lt:12+lt])[0]; o = 12 + lt
    for _ in range(nref):
        need(o + 4); ln = struct.unpack('<i', data[o:o+4])[0]; need(o + 8 + ln); o += 8 + ln
    del data[:o]
    while need(4):
        bs = struct.unpack('<i', data[:4])[0]
        need(4 + bs); yield bytes(data[4:4+bs]); del data[:4+bs]
n = 0; d = 0; first = None
ra, rb = records(sys.argv[1]), records(sys.argv[2])
while True:
    a = next(ra, None); b = next(rb, None)
    if a is None and b is None: break
    n += 1
    if a != b:
        d += 1
        if first is None: first = (n, a[:80] if a else None, b[:80] if b else None)
print(f"{n} records, {d} differ" + (f"; first at #{first[0]}" if first else ""))
