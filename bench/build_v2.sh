#!/bin/bash
# Build the x86-64-v2 (no compile-time AVX2) copy of the F015 raw-bit test
# binary for the SIMD axis of tests/f015.rs::f015_matrix.  In that build
# FOILEY_SIMD=0 forces the scalar kernels; the matrix compares it with the
# native (target-cpu=native, AVX2) build.  Output: target/x86-64-v2.
set -eu
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
# RUSTFLAGS replaces .cargo/config.toml's target-cpu=native
RUSTFLAGS="-C target-cpu=x86-64-v2" cargo test --release -p foiley --test f015 --no-run --target-dir target/x86-64-v2
