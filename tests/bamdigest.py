#!/usr/bin/env python3
"""bamdigest.py IN.bam: '<records> <md5>' of a BAM's records as bytes (block_size + record), header excluded.

Two BAMs with the same digest have the same records byte for byte, including the bin field and the encoding of
tags, which `samtools view` text would not show. Standard library only (a BGZF file is a multi-member gzip)."""
import gzip, hashlib, struct, sys


def digest(path):
    with gzip.open(path, 'rb') as f:
        if f.read(4) != b'BAM\1':
            raise SystemExit(f'{path}: not a BAM file')
        (l_text,) = struct.unpack('<i', f.read(4))
        f.read(l_text)
        (n_ref,) = struct.unpack('<i', f.read(4))
        for _ in range(n_ref):
            (l_name,) = struct.unpack('<i', f.read(4))
            f.read(l_name + 4)
        h, n = hashlib.md5(), 0
        while True:
            b = f.read(4)
            if not b:
                break
            (size,) = struct.unpack('<i', b)
            h.update(b)
            h.update(f.read(size))
            n += 1
    return n, h.hexdigest()


if __name__ == '__main__':
    for p in sys.argv[1:]:
        n, d = digest(p)
        print(f'{n} {d}')
