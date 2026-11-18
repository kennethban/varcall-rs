//! MD5 (RFC 1321), for the `.md5` files GATK / Picard write beside a BAM (`--create-output-bam-md5`,
//! `CREATE_MD5_FILE`): htsjdk's `Md5CalculatingOutputStream` writes the 32 lowercase hex digits, no newline.

const S: [u32; 64] = [
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11,
    16, 23, 4, 11, 16, 23, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
];

fn k() -> [u32; 64] {
    let mut k = [0u32; 64];
    for (i, x) in k.iter_mut().enumerate() {
        *x = ((i as f64 + 1.0).sin().abs() * 4294967296.0) as u32;
    }
    k
}

pub struct Md5 {
    state: [u32; 4],
    buf: Vec<u8>,
    len: u64,
    k: [u32; 64],
}

impl Default for Md5 {
    fn default() -> Self {
        Md5 { state: [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476], buf: Vec::with_capacity(64), len: 0, k: k() }
    }
}

impl Md5 {
    fn block(&mut self, b: &[u8]) {
        let mut m = [0u32; 16];
        for (i, w) in m.iter_mut().enumerate() {
            *w = u32::from_le_bytes(b[4 * i..4 * i + 4].try_into().unwrap());
        }
        let [mut a, mut bb, mut c, mut d] = self.state;
        for i in 0..64 {
            let (f, g) = match i / 16 {
                0 => ((bb & c) | (!bb & d), i),
                1 => ((d & bb) | (!d & c), (5 * i + 1) % 16),
                2 => (bb ^ c ^ d, (3 * i + 5) % 16),
                _ => (c ^ (bb | !d), (7 * i) % 16),
            };
            let f = f.wrapping_add(a).wrapping_add(self.k[i]).wrapping_add(m[g]);
            a = d;
            d = c;
            c = bb;
            bb = bb.wrapping_add(f.rotate_left(S[i]));
        }
        for (s, v) in self.state.iter_mut().zip([a, bb, c, d]) {
            *s = s.wrapping_add(v);
        }
    }

    pub fn update(&mut self, mut data: &[u8]) {
        self.len += data.len() as u64;
        if !self.buf.is_empty() {
            let take = (64 - self.buf.len()).min(data.len());
            self.buf.extend_from_slice(&data[..take]);
            data = &data[take..];
            if self.buf.len() < 64 {
                return;
            }
            let b = std::mem::take(&mut self.buf);
            self.block(&b);
            self.buf = b;
            self.buf.clear();
        }
        let mut chunks = data.chunks_exact(64);
        for b in &mut chunks {
            self.block(b);
        }
        self.buf.extend_from_slice(chunks.remainder());
    }

    pub fn hex(mut self) -> String {
        let bits = self.len.wrapping_mul(8);
        let mut pad = vec![0x80u8];
        while (self.buf.len() + pad.len()) % 64 != 56 {
            pad.push(0);
        }
        pad.extend_from_slice(&bits.to_le_bytes());
        let len = self.len;
        self.update(&pad);
        self.len = len;
        self.state.iter().flat_map(|s| s.to_le_bytes()).map(|b| format!("{b:02x}")).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::Md5;
    fn md5(parts: &[&[u8]]) -> String {
        let mut m = Md5::default();
        for p in parts {
            m.update(p);
        }
        m.hex()
    }
    #[test]
    fn rfc1321_vectors() {
        assert_eq!(md5(&[b""]), "d41d8cd98f00b204e9800998ecf8427e");
        assert_eq!(md5(&[b"abc"]), "900150983cd24fb0d6963f7d28e17f72");
        assert_eq!(md5(&[b"message ", b"digest"]), "f96b697d7cb7938d525a2f31aaf161d0");
        let s = b"12345678901234567890123456789012345678901234567890123456789012345678901234567890";
        assert_eq!(md5(&[&s[..3], &s[3..70], &s[70..]]), "57edf4a22be3c955ac49da2e2107b67a");
    }
}
