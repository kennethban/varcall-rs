#!/usr/bin/env python3
"""make.py IN.bam OUT.bam MODE: copy a BAM, setting some base quality bytes outside the SAM range 0-93 (which SAM
text cannot carry), in every 50th record:
  hi       bytes 200 / 130 / 200 at bases 40-42 (negative as Java bytes)
  noclip   the same, only in reads GATK's BaseRecalibrator does not clip (no S/H, a mate >= 400 bp away)
  mid      bytes 100 / 120 at bases 50-51 (94-127: GATK's tables have no cell for them)"""
import sys, struct, zlib
def bgzf_blocks(path):
    f = open(path, 'rb')
    while True:
        h = f.read(18)
        if len(h) < 18: return
        bsize = struct.unpack('<H', h[16:18])[0] + 1
        rest = f.read(bsize - 18)
        yield zlib.decompress(rest[:-8], -15)
def write_bgzf(out, data):
    for i in range(0, len(data), 0xff00):
        chunk = data[i:i + 0xff00]
        c = zlib.compressobj(6, zlib.DEFLATED, -15); z = c.compress(chunk) + c.flush()
        out.write(struct.pack('<4BI2BH2BH', 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2) + struct.pack('<H', len(z) + 25) + z
                  + struct.pack('<II', zlib.crc32(chunk) & 0xffffffff, len(chunk)))
src, dst, mode = sys.argv[1], sys.argv[2], sys.argv[3]
data = b''.join(bgzf_blocks(src))
lt = struct.unpack('<i', data[4:8])[0]; o = 8 + lt; nref = struct.unpack('<i', data[o:o + 4])[0]; o += 4
for _ in range(nref):
    ln = struct.unpack('<i', data[o:o + 4])[0]; o += 8 + ln
out = bytearray(data[:o]); n = 0; edited = 0
while o < len(data):
    bs = struct.unpack('<i', data[o:o + 4])[0]; rec = bytearray(data[o + 4:o + 4 + bs]); o += 4 + bs
    l_rn, n_cig, l_seq = rec[8], struct.unpack('<H', rec[12:14])[0], struct.unpack('<i', rec[16:20])[0]
    q0 = 32 + l_rn + 4 * n_cig + (l_seq + 1) // 2
    cig = [struct.unpack('<I', rec[32 + l_rn + 4 * k:36 + l_rn + 4 * k])[0] & 0xf for k in range(n_cig)]
    flag, tlen = struct.unpack('<H', rec[14:16])[0], struct.unpack('<i', rec[28:32])[0]
    clipped = any(c in (4, 5) for c in cig) or not (flag & 1) or abs(tlen) < 400 or (flag & 8)
    if n % 50 == 0 and l_seq > 60 and rec[q0] != 0xff and not (mode == 'noclip' and clipped):
        if mode in ('hi', 'noclip'):
            rec[q0 + 40] = 200; rec[q0 + 41] = 130; rec[q0 + 42] = 200
        else:
            rec[q0 + 50] = 100; rec[q0 + 51] = 120
        edited += 1
    n += 1
    out += struct.pack('<i', bs) + rec
with open(dst, 'wb') as f:
    write_bgzf(f, bytes(out))
    f.write(bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000'))
print(f'{n} records, {edited} edited')
