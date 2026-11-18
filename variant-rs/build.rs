// hc-cudahmm's build script sets the rpath only for its own artifacts (cargo does not propagate
// `rustc-link-arg` from dependencies), so the binaries of this crate set it again.
// CRAM input: csrc/cram_shim.c against htslib from the conda env (HTSLIB_PREFIX overrides).
fn main() {
    let m = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let dir = std::env::var("CUDA_PAIRHMM_DIR").unwrap_or_else(|_| format!("{m}/cuda-pairhmm/build"));
    if let Ok(dir) = std::fs::canonicalize(&dir) {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}", dir.display());
    }
    let hts = std::env::var("HTSLIB_PREFIX").or_else(|_| std::env::var("CONDA_PREFIX")).expect("set HTSLIB_PREFIX or activate the conda env with htslib");
    cc::Build::new().file("csrc/cram_shim.c").include(format!("{hts}/include")).opt_level(2).compile("cram_shim");
    println!("cargo:rustc-link-search=native={hts}/lib");
    println!("cargo:rustc-link-lib=dylib=hts");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{hts}/lib");
    println!("cargo:rerun-if-changed=csrc/cram_shim.c");
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
    println!("cargo:rerun-if-env-changed=HTSLIB_PREFIX");
    println!("cargo:rerun-if-env-changed=CUDA_PAIRHMM_DIR");
}
