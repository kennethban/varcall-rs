//! FASTQ input as bwa-mem2 reads it (kseq.h `kseq_read`, bwa.cpp `bseq_read_orig`): batches of at least `chunk`
//! bases with an even number of reads (paired: R1, R2 interleaved), names cut at the first whitespace and a
//! trailing "/1" or "/2" removed, bases kept as read (encoded to 0..4 later).

use std::io::{BufRead, BufReader, Read};

pub struct Seq {
    pub name: Vec<u8>,
    pub comment: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

pub struct Fastq {
    r: Box<dyn BufRead + Send>,
    /// kseq's last_char: the next record's '@' / '>' already read (0: not yet)
    last_char: u8,
}

fn open_reader(path: &str) -> Box<dyn BufRead + Send> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut b = BufReader::with_capacity(1 << 20, f);
    let gz = b.fill_buf().map(|h| h.len() >= 2 && h[0] == 0x1f && h[1] == 0x8b).unwrap_or(false);
    if gz {
        // decompressed on a thread of its own, handed over in 4 MB blocks: the read stage only parses
        Box::new(Threaded::spawn(GzMembers { dec: Some(flate2::bufread::GzDecoder::new(b)) }))
    } else {
        Box::new(b)
    }
}

/// Concatenated gzip members, as zlib's gzread reads them: after a member, bytes that do not start another gzip
/// header (zero padding, a stray newline) end the input instead of being an error (gz_look).
struct GzMembers {
    dec: Option<flate2::bufread::GzDecoder<BufReader<std::fs::File>>>,
}

impl Read for GzMembers {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        loop {
            let Some(d) = self.dec.as_mut() else { return Ok(0) };
            let n = d.read(out)?;
            if n > 0 || out.is_empty() {
                return Ok(n);
            }
            // end of a member: another one, or the end of the input
            let mut inner = self.dec.take().unwrap().into_inner();
            let h = inner.fill_buf()?;
            if h.first() == Some(&0x1f) && h.get(1).is_none_or(|&c| c == 0x8b) {
                self.dec = Some(flate2::bufread::GzDecoder::new(inner));
            }
        }
    }
}

/// A reader running on its own thread: blocks of its output arrive through a bounded channel
struct Threaded {
    rx: std::sync::mpsc::Receiver<std::io::Result<Vec<u8>>>,
    buf: Vec<u8>,
    pos: usize,
    done: bool,
}

impl Threaded {
    fn spawn<R: Read + Send + 'static>(mut r: R) -> Threaded {
        let (tx, rx) = std::sync::mpsc::sync_channel::<std::io::Result<Vec<u8>>>(4);
        std::thread::Builder::new()
            .name("gunzip".into())
            .spawn(move || loop {
                let mut v = vec![0u8; 4 << 20];
                let mut n = 0;
                // fill the block (a read returns at most one deflate stretch)
                let res = loop {
                    match r.read(&mut v[n..]) {
                        Ok(0) => break Ok(()),
                        Ok(k) => {
                            n += k;
                            if n == v.len() {
                                break Ok(());
                            }
                        }
                        Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                        Err(e) => break Err(e),
                    }
                };
                v.truncate(n);
                let end = n == 0;
                let msg = match res {
                    Ok(()) => Ok(v),
                    Err(e) => Err(e),
                };
                let err = msg.is_err();
                if tx.send(msg).is_err() || end || err {
                    return;
                }
            })
            .expect("spawning the gunzip thread");
        Threaded { rx, buf: Vec::new(), pos: 0, done: false }
    }
}

impl Read for Threaded {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        let b = self.fill_buf()?;
        let n = b.len().min(out.len());
        out[..n].copy_from_slice(&b[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for Threaded {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        while self.pos == self.buf.len() && !self.done {
            match self.rx.recv() {
                Ok(Ok(v)) if v.is_empty() => self.done = true,
                Ok(Ok(v)) => {
                    self.buf = v;
                    self.pos = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => self.done = true,
            }
        }
        Ok(&self.buf[self.pos..])
    }
    fn consume(&mut self, n: usize) {
        self.pos += n;
    }
}

impl Fastq {
    pub fn open(path: &str) -> Fastq {
        Fastq { r: open_reader(path), last_char: 0 }
    }

    /// ks_getc
    #[inline]
    fn getc(&mut self) -> Option<u8> {
        let b = self.r.fill_buf().expect("reading FASTQ");
        let c = *b.first()?;
        self.r.consume(1);
        Some(c)
    }

    /// ks_getuntil2 up to the first byte `stop` accepts, appended to `out`; that byte is consumed. None at the end of
    /// the input with nothing read (-1), else the stopping byte (None if the input ended). `line`: KS_SEP_LINE, which
    /// drops a final '\r' when the whole string is longer than 1.
    fn get_until(&mut self, out: &mut Vec<u8>, line: bool) -> Option<Option<u8>> {
        let mut gotany = false;
        let mut stop = None;
        loop {
            let b = self.r.fill_buf().expect("reading FASTQ");
            if b.is_empty() {
                break;
            }
            gotany = true;
            let at = if line { memchr::memchr(b'\n', b) } else { b.iter().position(|&c| is_space(c)) };
            match at {
                Some(i) => {
                    out.extend_from_slice(&b[..i]);
                    stop = Some(b[i]);
                    self.r.consume(i + 1);
                    break;
                }
                None => {
                    let n = b.len();
                    out.extend_from_slice(b);
                    self.r.consume(n);
                }
            }
        }
        if !gotany {
            return None;
        }
        if line && out.len() > 1 && out.last() == Some(&b'\r') {
            out.pop();
        }
        Some(stop)
    }

    /// kseq_read (kseq.h), byte for byte: Ok(None) at the end of the input (-1), Err(()) for a FASTQ record without a
    /// quality line or with one of a different length than its sequence (-2; the bytes read are consumed).
    pub fn next(&mut self) -> Result<Option<Seq>, ()> {
        if self.last_char == 0 {
            // jump to the next header: any byte up to a '>' or '@' is skipped
            loop {
                match self.getc() {
                    None => return Ok(None),
                    Some(c @ (b'>' | b'@')) => {
                        self.last_char = c;
                        break;
                    }
                    Some(_) => {}
                }
            }
        }
        let mut name = Vec::new();
        let Some(c) = self.get_until(&mut name, false) else { return Ok(None) };
        let mut comment = Vec::new();
        if c != Some(b'\n') {
            self.get_until(&mut comment, true);
        }
        // trim_readno
        if name.len() > 2 && name[name.len() - 2] == b'/' && name[name.len() - 1].is_ascii_digit() {
            name.truncate(name.len() - 2);
        }
        let mut seq = Vec::new();
        let mut c;
        loop {
            c = self.getc();
            match c {
                None | Some(b'>' | b'+' | b'@') => break,
                Some(b'\n') => continue, // empty lines
                Some(x) => {
                    seq.push(x);
                    self.get_until(&mut seq, true);
                }
            }
        }
        if let Some(h @ (b'>' | b'@')) = c {
            self.last_char = h;
        }
        if c != Some(b'+') {
            // FASTA (or a FASTQ record cut before its '+')
            return Ok(Some(Seq { name, comment, seq, qual: Vec::new() }));
        }
        // the rest of the '+' line
        loop {
            match self.getc() {
                None => return Err(()),
                Some(b'\n') => break,
                Some(_) => {}
            }
        }
        // quality lines: one, then more while shorter than the sequence
        let mut qual = Vec::new();
        while self.get_until(&mut qual, true).is_some() && qual.len() < seq.len() {}
        self.last_char = 0;
        if qual.len() != seq.len() {
            return Err(());
        }
        Ok(Some(Seq { name, comment, seq, qual }))
    }
}

/// C isspace in the "C" locale
#[inline]
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// bseq_read_orig: the next batch (paired: R1, R2, R1, R2 ...), at least `chunk` bases and an even count. A
/// malformed record ends the batch and is dropped (`while (kseq_read(ks) >= 0)`), with its mate; reading goes on
/// after it in the next batch.
pub fn read_batch(r1: &mut Fastq, r2: Option<&mut Fastq>, chunk: i64) -> Vec<Seq> {
    let mut out = Vec::new();
    let mut size = 0i64;
    let mut r2 = r2;
    while let Ok(Some(a)) = r1.next() {
        let b = match r2.as_mut() {
            Some(r) => match r.next() {
                Ok(Some(b)) => Some(b),
                _ => {
                    eprintln!("[W::bseq_read_orig] the 2nd file has fewer sequences.");
                    break;
                }
            },
            None => None,
        };
        size += a.seq.len() as i64;
        out.push(a);
        if let Some(b) = b {
            size += b.seq.len() as i64;
            out.push(b);
        }
        if size >= chunk && out.len() % 2 == 0 {
            break;
        }
    }
    if size == 0 {
        // test if the 2nd file is finished (reads, and drops, one of its records)
        if let Some(r) = r2.as_mut() {
            if let Ok(Some(_)) = r.next() {
                eprintln!("[W::bseq_read_orig] the 1st file has fewer sequences.");
            }
        }
    }
    out
}

/// nst_nt4_table
pub fn nt4(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

#[allow(dead_code)]
fn _assert_read<R: Read>(_: R) {}
