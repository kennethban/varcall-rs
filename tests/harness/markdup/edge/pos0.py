#!/usr/bin/env python3
"""pos0.py IN.bam OUT.bam: copy a coordinate-sorted BAM, adding before its first record an unmapped copy of it that
keeps the contig but has POS 0 (SAM allows it; htslib's SAM parser drops such an RNAME, so this writes BAM bytes).
Review finding 32: Picard / GATK copy such records through; variant-rs dropped them."""
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
src, dst = sys.argv[1], sys.argv[2]
data = b''.join(bgzf_blocks(src))
lt = struct.unpack('<i', data[4:8])[0]; o = 8 + lt; nref = struct.unpack('<i', data[o:o + 4])[0]; o += 4
for _ in range(nref):
    ln = struct.unpack('<i', data[o:o + 4])[0]; o += 8 + ln
out = bytearray(data[:o])
bs = struct.unpack('<i', data[o:o + 4])[0]; rec = data[o + 4:o + 4 + bs]
tid = struct.unpack('<i', rec[0:4])[0]
l_rn, n_cig, l_seq = rec[8], struct.unpack('<H', rec[12:14])[0], struct.unpack('<i', rec[16:20])[0]
name = rec[32:32 + l_rn - 1] + b'_pos0\0'
rest = rec[32 + l_rn + 4 * n_cig:32 + l_rn + 4 * n_cig + (l_seq + 1) // 2 + l_seq] + rec[32 + l_rn + 4 * n_cig + (l_seq + 1) // 2 + l_seq:]
# refID, pos -1, bin 4680 (reg2bin(-1, 0)) | mapq 0 | l_read_name, flag 4 | n_cigar 0, l_seq, next refID -1, next pos -1, tlen 0
new = struct.pack('<iiIIiiii', tid, -1, (4680 << 16) | len(name), (4 << 16), l_seq, -1, -1, 0) + name + rest
out += struct.pack('<i', len(new)) + new + data[o:]
with open(dst, 'wb') as f:
    write_bgzf(f, bytes(out))
    f.write(bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000'))
