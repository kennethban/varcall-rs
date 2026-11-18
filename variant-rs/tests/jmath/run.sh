#!/bin/bash
# Compare variant-rs/crates/assembly/src/jmath.rs (+ fdlibm.rs) with the running JVM's Math.exp/log/log10/pow/log1p/expm1, bit for bit.
#   tests/jmath/run.sh [millions of inputs per function, default 10] [seed]
# Needs `java` (17) and `rustc` on PATH.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
out=${TMPDIR:-/tmp}/jmath_check.$$
rustc -O -C target-cpu=native ${CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER:+-C linker=$CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER} --edition 2021 -o "$out" "$here/check.rs"
trap 'rm -f "$out"' EXIT
java -version 2>&1 | head -1
java "$here/JMathGen.java" "${1:-10}" "${2:-20261116}" | "$out"
