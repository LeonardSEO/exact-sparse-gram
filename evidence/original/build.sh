#!/bin/sh
set -eu
cd "$(dirname "$0")"
CFLAGS='-O3 -std=c11 -fno-fast-math -ffp-contract=off -Wall -Wextra'
SOURCES='src/joint_kernel.c src/communication_skinny.c src/adapters.c src/algebra_hankel.c src/numerical_guard.c src/algebra_sparse.c src/numerical_sparse.c src/communication_guard.c src/communication_guard_fast.c src/numerical_sparse_fast.c'
clang $CFLAGS src/benchmark.c $SOURCES -framework Accelerate -o benchmark
LATTICE='src/lattice_exact.c src/lattice_substitute.c src/lattice_sequential.c'
clang $CFLAGS src/lattice_benchmark.c $SOURCES $LATTICE -framework Accelerate -o lattice_benchmark
clang $CFLAGS -shared -fPIC $SOURCES $LATTICE -o libgram.dylib
