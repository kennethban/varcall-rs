// Links libcudapairhmm.so (make -C variant-rs/cuda-pairhmm). Without it (no CUDA toolkit), the crate builds with
// cfg no_cuda: init() reports no GPU and variant-rs runs the CPU PairHMM, which gives the same results.
fn main() {
    println!("cargo::rustc-check-cfg=cfg(no_cuda)");
    let dir = std::env::var("CUDA_PAIRHMM_DIR").unwrap_or_else(|_| {
        let m = std::env::var("CARGO_MANIFEST_DIR").unwrap();
        format!("{m}/../../cuda-pairhmm/build")
    });
    println!("cargo:rerun-if-env-changed=CUDA_PAIRHMM_DIR");
    println!("cargo:rerun-if-changed={dir}/libcudapairhmm.so");
    let Ok(dir) = std::fs::canonicalize(&dir).map_err(|_| ()).and_then(|d| if d.join("libcudapairhmm.so").exists() { Ok(d) } else { Err(()) }) else {
        println!("cargo:warning=libcudapairhmm.so not built: CPU PairHMM only (make -C variant-rs/cuda-pairhmm for the GPU one)");
        println!("cargo:rustc-cfg=no_cuda");
        return;
    };
    println!("cargo:rustc-link-search=native={}", dir.display());
    println!("cargo:rustc-link-lib=dylib=cudapairhmm");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", dir.display());
    // libcudapairhmm.so needs libstdc++ from the conda toolchain env
    if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
        println!("cargo:rustc-link-search=native={prefix}/lib");
        println!("cargo:rustc-link-arg=-Wl,-rpath,{prefix}/lib");
        println!("cargo:rustc-link-lib=dylib=stdc++");
    }
}
