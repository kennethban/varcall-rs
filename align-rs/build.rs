// Links cuda/build/libalngpu.so (make -C cuda) when it exists; without it align-rs builds CPU-only.
fn main() {
    let m = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    // hts/alnhts.c against htslib (SAM -> BAM records, the sorted BAM writer and its index): HTSLIB_PREFIX, else
    // the active conda env
    let prefix = std::env::var("HTSLIB_PREFIX").or_else(|_| std::env::var("CONDA_PREFIX")).expect("set HTSLIB_PREFIX or activate the conda env with htslib");
    println!("cargo:rerun-if-env-changed=HTSLIB_PREFIX");
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
    println!("cargo:rerun-if-changed={m}/hts/alnhts.c");
    cc::Build::new().file(format!("{m}/hts/alnhts.c")).include(format!("{prefix}/include")).opt_level(2).compile("alnhts");
    println!("cargo:rustc-link-search=native={prefix}/lib");
    println!("cargo:rustc-link-lib=dylib=hts");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{prefix}/lib");
    let dir = format!("{m}/cuda/build");
    println!("cargo:rerun-if-changed={dir}/libalngpu.so");
    if std::path::Path::new(&format!("{dir}/libalngpu.so")).exists() {
        println!("cargo:rustc-link-search=native={dir}");
        println!("cargo:rustc-link-lib=dylib=alngpu");
        println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
        if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
            println!("cargo:rustc-link-search=native={prefix}/lib");
            println!("cargo:rustc-link-arg=-Wl,-rpath,{prefix}/lib");
            println!("cargo:rustc-link-lib=dylib=stdc++");
        }
        println!("cargo:rustc-cfg=gpu");
    }
    println!("cargo::rustc-check-cfg=cfg(gpu)");
}
