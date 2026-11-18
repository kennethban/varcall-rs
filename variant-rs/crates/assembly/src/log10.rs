//! Java's Math.log10 (HotSpot x86_64 intrinsic) differs from glibc's log10 in the last
//! ulp for ~1.6% of integers. K-best path scores are sums of log10 differences of edge
//! multiplicities, so we use a table of values produced by Java itself for 0..2^17.

static TABLE_BYTES: &[u8] = include_bytes!("log10_java.bin");

#[inline]
pub fn java_log10(x: i32) -> f64 {
    let i = x as usize;
    if x >= 0 && i * 8 + 8 <= TABLE_BYTES.len() {
        let b: [u8; 8] = TABLE_BYTES[i * 8..i * 8 + 8].try_into().unwrap();
        f64::from_le_bytes(b)
    } else {
        crate::jmath::log10(x as f64)
    }
}
