// hc-cudahmm's build script sets the rpath only for its own artifacts (cargo does not propagate
// `rustc-link-arg` from dependencies), so the binaries of this crate set it again.
fn main() {
    let m = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let dir = std::env::var("CUDA_PAIRHMM_DIR").unwrap_or_else(|_| format!("{m}/../../cuda-pairhmm/build"));
    if let Ok(dir) = std::fs::canonicalize(&dir) {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}", dir.display());
    }
    if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{prefix}/lib");
    }
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
    println!("cargo:rerun-if-env-changed=CUDA_PAIRHMM_DIR");
}
