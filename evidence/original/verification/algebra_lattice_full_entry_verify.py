#!/usr/bin/env python3
"""Independent full-entry exact oracle for held-out final lattice cases.

Default: print the reproduction plan without compiling or running anything.
After the benchmark owner explicitly releases its execution lock, run:
  python3 gram_research/verification/algebra_lattice_full_entry_verify.py --run
This is untimed verification: the C executable has no timer or benchmark loop.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

BASE = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
CHECKER = HERE / "algebra_lattice_full_entry.c"
KERNEL = BASE / "src/lattice_sequential.c"
TRACKED = [CHECKER, Path(__file__), KERNEL,
           BASE / "src/lattice_benchmark.c", BASE / "src/benchmark.c"]
BINARY = HERE / "algebra_lattice_full_entry"
EVENTS = HERE / "algebra_lattice_full_entry_errors.tsv"
RECEIPT = HERE / "algebra_lattice_full_entry_results.json"
RAW = HERE / "algebra_lattice_full_entry_output.jsonl"
TAU = Fraction.from_float(2e-12)


def hashes():
    return {str(p.relative_to(BASE)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in TRACKED}


def dyadic(exponent):
    return Fraction(1 << exponent) if exponent >= 0 else Fraction(1, 1 << -exponent)


def check_exact_tolerances():
    count = failures = new_errors = 0
    by_method = {}
    with EVENTS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            method = row["method"]
            got = Fraction.from_float(float.fromhex(row["got_hex"]))
            integer = int(row["integer_target"])
            exponent = int(row["scale_exp"])
            exact = integer * dyadic(2 * exponent)
            error = abs(got - exact)
            norm_product_square = (int(row["integer_norm_i"]) *
                                   int(row["integer_norm_j"]) * dyadic(4 * exponent))
            passed = error * error <= TAU * TAU * norm_product_square
            method_state = by_method.setdefault(method,
                {"nonzero_errors": 0, "preset_tolerance_failures": 0,
                 "maximum_normalized_squared_error_exact": Fraction(0)})
            method_state["nonzero_errors"] += 1
            method_state["preset_tolerance_failures"] += not passed
            if norm_product_square:
                ratio_squared = error * error / norm_product_square
                method_state["maximum_normalized_squared_error_exact"] = max(
                    method_state["maximum_normalized_squared_error_exact"], ratio_squared)
            elif error:
                method_state["zero_norm_error"] = True
            count += 1
            failures += not passed
            new_errors += method == "lattice_sequential" and error != 0
    for state in by_method.values():
        state["maximum_normalized_squared_error_exact"] = str(
            state["maximum_normalized_squared_error_exact"])
    return {"lossless_nonzero_error_events": count,
            "preset_tolerance_failures_exact": failures,
            "new_kernel_nonzero_errors_exact": new_errors,
            "nonzero_error_events_by_method": by_method}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true",
        help="Compile and run only after the benchmark owner releases its lock")
    args = parser.parse_args()
    command = ["clang", "-O3", "-std=c11", "-fno-fast-math", "-ffp-contract=off",
               "-Wall", "-Wextra", "-Wconversion", "-Werror",
               str(CHECKER), str(KERNEL), "-framework", "Accelerate", "-o", str(BINARY)]
    if not args.run:
        print(json.dumps({"state": "prepared_only_no_compile_or_run",
                          "compile_command": command,
                          "run_command": [str(BINARY), str(EVENTS)],
                          "receipt": str(RECEIPT)}, indent=2))
        return
    if sys.platform != "darwin":
        raise SystemExit("This receipt compares Apple Accelerate and requires macOS.")
    before = hashes()
    compilation = subprocess.run(command, capture_output=True, text=True)
    receipt = {"artifact": "untimed independent full-entry integer Gram oracle",
        "compile_command": command, "compile_exit": compilation.returncode,
        "compile_stdout": compilation.stdout, "compile_stderr": compilation.stderr,
        "source_sha256_before": before,
        "platform": platform.platform(), "machine": platform.machine(),
        "python": sys.version, "preset_tolerance_binary64_exact": str(TAU),
        "oracle": "Every lower entry: direct int32 products accumulated in int64, absolute bound <=2^46; target converted exactly by power-of-two scaling",
        "data": "Integer-first xorshift64 generator mirrors lattice_benchmark.c generation order; no candidate or benchmark generator is called",
        "compiler_version": subprocess.run(["clang", "--version"],capture_output=True,text=True).stdout,
        "timing": "No clocks, elapsed times, repetitions, or benchmark measurements are used"}
    if compilation.returncode:
        receipt["passed"] = False
        RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n")
        raise SystemExit(compilation.stderr)
    env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS="1")
    execution = subprocess.run([str(BINARY), str(EVENTS)],capture_output=True,text=True,env=env)
    RAW.write_text(execution.stdout)
    receipt.update({"run_command": [str(BINARY), str(EVENTS)],
                    "run_exit": execution.returncode, "run_stderr": execution.stderr,
                    "cases_and_layouts": [json.loads(line) for line in execution.stdout.splitlines() if line.strip()],
                    "binary_sha256": hashlib.sha256(BINARY.read_bytes()).hexdigest(),
                    "events_sha256": hashlib.sha256(EVENTS.read_bytes()).hexdigest() if EVENTS.exists() else None,
                    "source_sha256_after": hashes()})
    if EVENTS.exists(): receipt.update(check_exact_tolerances())
    complete = len(receipt["cases_and_layouts"]) == 8
    receipt["passed"] = (compilation.returncode == 0 and execution.returncode == 0 and
                         complete and before == receipt["source_sha256_after"] and
                         receipt.get("preset_tolerance_failures_exact") == 0 and
                         receipt.get("new_kernel_nonzero_errors_exact") == 0)
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    if not receipt["passed"]: raise SystemExit(1)


if __name__ == "__main__": main()
