//! CRAM input through htslib (`../csrc/cram_shim.c`, linked against the conda env's libhts).
//!
//! Records are decoded by htslib, the library that also converts CRAM to BAM (`samtools view -b`), and
//! serialized as BAM record bytes, so the front end sees exactly the records of the converted BAM.
//! A pool of open handles (each with its own index, iterator and reference cache) is reused across contigs.

use std::ffi::{c_char, c_int, CStr, CString};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{sync_channel, Receiver};
use std::sync::{Arc, Mutex};

#[repr(C)]
struct Crh {
    _p: [u8; 0],
}

extern "C" {
    fn crh_open(path: *const c_char, fasta: *const c_char, share: *mut Crh) -> *mut Crh;
    fn crh_close(h: *mut Crh);
    fn crh_header_text(h: *mut Crh, len: *mut usize) -> *const c_char;
    fn crh_nref(h: *mut Crh) -> c_int;
    fn crh_ref_name(h: *mut Crh, tid: c_int) -> *const c_char;
    fn crh_ref_len(h: *mut Crh, tid: c_int) -> i64;
    fn crh_query(h: *mut Crh, tid: c_int, iv: *const i64, n: c_int) -> c_int;
    fn crh_query_unplaced(h: *mut Crh) -> c_int;
    fn crh_next(h: *mut Crh, out: *mut *const u8) -> i64;
}

struct Handle(*mut Crh);
// a handle is used by one thread at a time (taken from and returned to the pool)
unsafe impl Send for Handle {}

pub struct CramPool {
    path: CString,
    fasta: CString,
    /// the handle opened first (header)
    first: Handle,
    free: Mutex<Vec<Handle>>,
}
unsafe impl Sync for CramPool {}
unsafe impl Send for CramPool {}

impl CramPool {
    pub fn open(path: &str, fasta: &str) -> CramPool {
        let (p, f) = (CString::new(path).unwrap(), CString::new(fasta).unwrap());
        let h = unsafe { crh_open(p.as_ptr(), f.as_ptr(), std::ptr::null_mut()) };
        assert!(!h.is_null(), "{path}: cannot open as CRAM with its index (.crai) and reference {fasta}");
        CramPool { path: p, fasta: f, first: Handle(h), free: Mutex::new(Vec::new()) }
    }

    fn take(&self) -> Handle {
        if let Some(h) = self.free.lock().unwrap().pop() {
            return h;
        }
        let h = unsafe { crh_open(self.path.as_ptr(), self.fasta.as_ptr(), self.first.0) };
        assert!(!h.is_null(), "cannot reopen {:?}", self.path);
        Handle(h)
    }

    /// (header text, reference names, reference lengths)
    pub fn header(&self) -> (String, Vec<String>, Vec<i64>) {
        let h = self.first.0;
        unsafe {
            let mut len = 0usize;
            let t = crh_header_text(h, &mut len);
            let text = String::from_utf8_lossy(std::slice::from_raw_parts(t as *const u8, len)).trim_end_matches('\0').to_string();
            let n = crh_nref(h);
            let names = (0..n).map(|i| CStr::from_ptr(crh_ref_name(h, i)).to_string_lossy().into_owned()).collect();
            let lens = (0..n).map(|i| crh_ref_len(h, i)).collect();
            (text, names, lens)
        }
    }

    /// Streams the records of `tid` overlapping `intervals` (1-based inclusive, sorted, disjoint), in file
    /// order, decoded on a reader thread.
    pub fn stream(self: &Arc<Self>, tid: usize, intervals: &[(i64, i64)], decode_ns: Arc<AtomicU64>) -> RecordStream {
        self.stream_q(Some(tid), intervals, decode_ns)
    }

    /// Streams the unplaced unmapped records (RNAME '*'), in file order
    pub fn stream_unplaced(self: &Arc<Self>, decode_ns: Arc<AtomicU64>) -> RecordStream {
        self.stream_q(None, &[], decode_ns)
    }

    fn stream_q(self: &Arc<Self>, tid: Option<usize>, intervals: &[(i64, i64)], decode_ns: Arc<AtomicU64>) -> RecordStream {
        let (tx, rx) = sync_channel::<Batch>(4);
        let pool = self.clone();
        let iv: Vec<i64> = intervals.iter().flat_map(|&(s, e)| [s, e]).collect();
        std::thread::Builder::new()
            .name("cram-read".into())
            .spawn(move || {
                let h = pool.take();
                let t0 = std::time::Instant::now();
                let rc = match tid {
                    Some(t) => unsafe { crh_query(h.0, t as c_int, iv.as_ptr(), (iv.len() / 2) as c_int) },
                    None => unsafe { crh_query_unplaced(h.0) },
                };
                assert_eq!(rc, 0, "CRAM query failed (tid {tid:?})");
                let mut b = Batch::default();
                let mut ptr: *const u8 = std::ptr::null();
                let mut ok = true;
                loop {
                    let n = unsafe { crh_next(h.0, &mut ptr) };
                    if n == -1 {
                        break;
                    }
                    assert!(n >= 0, "CRAM decoding failed (tid {tid:?}, code {n})");
                    b.data.extend_from_slice(unsafe { std::slice::from_raw_parts(ptr, n as usize) });
                    b.ends.push(b.data.len() as u32);
                    if b.data.len() >= (1 << 20) {
                        if tx.send(std::mem::take(&mut b)).is_err() {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    b.end = true; // a stream that closes without this marker failed
                    let _ = tx.send(b);
                }
                decode_ns.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
                pool.free.lock().unwrap().push(h);
            })
            .unwrap();
        RecordStream { rx, cur: Batch::default(), i: 0 }
    }
}

impl Drop for CramPool {
    fn drop(&mut self) {
        for h in self.free.get_mut().unwrap().drain(..) {
            unsafe { crh_close(h.0) };
        }
        unsafe { crh_close(self.first.0) };
    }
}

#[derive(Default)]
struct Batch {
    data: Vec<u8>,
    ends: Vec<u32>,
    end: bool,
}

pub struct RecordStream {
    rx: Receiver<Batch>,
    cur: Batch,
    i: usize,
}

impl RecordStream {
    fn done(&self) -> bool {
        self.cur.end
    }
}

impl RecordStream {
    /// Next record (the bytes after `block_size`) into `out`. False at the end.
    pub fn next_record(&mut self, out: &mut Vec<u8>) -> bool {
        while self.i >= self.cur.ends.len() {
            if self.done() {
                return false;
            }
            match self.rx.recv() {
                Ok(b) => {
                    self.cur = b;
                    self.i = 0;
                }
                Err(_) => panic!("CRAM reader thread failed"),
            }
        }
        let s = if self.i == 0 { 0 } else { self.cur.ends[self.i - 1] as usize };
        let e = self.cur.ends[self.i] as usize;
        out.clear();
        out.extend_from_slice(&self.cur.data[s..e]);
        self.i += 1;
        true
    }
}

pub fn crai_path(cram: &str) -> String {
    [format!("{cram}.crai"), cram.strip_suffix(".cram").map(|s| format!("{s}.crai")).unwrap_or_default()]
        .into_iter()
        .find(|p| std::path::Path::new(p).is_file())
        .unwrap_or_else(|| panic!("no .crai index for {cram}"))
}

/// Per-16 kb-window cost from the CRAI (compressed slice bytes spread evenly over each slice's span),
/// the CRAM counterpart of `plan::window_costs` on a BAI.
pub fn crai_window_costs(crai: &str, lens: &[i64]) -> Vec<Vec<f64>> {
    use std::io::Read;
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(crai).unwrap_or_else(|e| panic!("{crai}: {e}")))
        .read_to_string(&mut text)
        .unwrap_or_else(|e| panic!("{crai}: {e}"));
    let w = crate::plan::WIN;
    let mut costs: Vec<Vec<f64>> = lens.iter().map(|&l| vec![0.0; ((l + w - 1) / w) as usize]).collect();
    for line in text.lines() {
        let f: Vec<i64> = line.split('\t').filter_map(|x| x.trim().parse().ok()).collect();
        if f.len() < 6 || f[0] < 0 || f[0] as usize >= costs.len() || f[2] <= 0 {
            continue; // unmapped or multi-reference slices
        }
        let v = &mut costs[f[0] as usize];
        if v.is_empty() {
            continue;
        }
        let (s, e) = (f[1].max(1), f[1] + f[2] - 1); // 1-based inclusive span
        let (ws, we) = (((s - 1) / w) as usize, (((e - 1) / w) as usize).min(v.len() - 1));
        let ws = ws.min(we);
        let per = f[5] as f64 / (we - ws + 1) as f64;
        for x in &mut v[ws..=we] {
            *x += per;
        }
    }
    costs
}
