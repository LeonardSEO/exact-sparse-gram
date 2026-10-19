#!/usr/bin/env python3
"""Untimed exact checks for the sequential-factor lattice Gram kernel.

No benchmark, timer, BLAS, dense conversion, or approximate reference is used.
The direct Gram reference uses Python integers and exact Fraction scaling.
"""
import ctypes as C
import ctypes.util
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / "src/lattice_sequential.c"
ORIGINAL = BASE / "src/lattice_exact.c"
ORIGINAL_SHA256 = "3402d7c6531fdb17178493df4205b54aaa580fcc6f6d99aa3a66cffa3e3da540"
SUBSTITUTE = BASE / "src/lattice_substitute.c"
SUBSTITUTE_SHA256 = "d1400c003886f0e57102e75bb7e442b7329cbad26d6d796b844b9bf877f1ba95"
PTR = C.POINTER(C.c_double)
SENTINEL = math.nan
LIMIT = 1 << 53


def build(scalar):
    suffix = "scalar" if scalar else "native"
    library = BASE / f"verification/lattice_sequential_{suffix}.dylib"
    cmd = ["clang", "-O3", "-std=c11", "-fno-fast-math", "-ffp-contract=off",
           "-Wall", "-Wextra", "-Wconversion", "-Werror", "-dynamiclib",
           str(SOURCE), "-o", str(library)]
    if scalar:
        cmd.insert(1, "-DLATTICE_SEQUENTIAL_SCALAR")
    if sys.platform != "darwin":
        cmd[cmd.index("-dynamiclib")] = "-shared"
        cmd.extend(["-fPIC", "-lm"])
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    function = C.CDLL(str(library)).lattice_sequential_gram
    function.argtypes = [PTR, C.c_size_t, C.c_size_t, C.c_size_t, PTR,
                         C.c_size_t, C.c_int, C.c_size_t]
    function.restype = C.c_int
    return function, cmd


def exact_scale(integer, exponent):
    return Fraction(integer << exponent) if exponent >= 0 else \
           Fraction(integer, 1 << -exponent)


def cap(n):
    return math.isqrt(LIMIT // n)


def shifted(rng, m, n, q, bound):
    rows = [[rng.randint(-bound, bound) for _ in range(n)]]
    for _ in range(1, m):
        row = rows[-1][1:] + [rng.randint(-bound, bound)]
        for k in rng.sample(range(n - 1), min(q, n - 1)):
            row[k] = rng.randint(-bound, bound)
        rows.append(row)
    return rows


def check_case(functions, name, rows, exponent, max_q):
    m, n = len(rows), len(rows[0])
    bound = cap(n)
    assert n * bound * bound <= LIMIT < n * (bound + 1) * (bound + 1)
    assert all(abs(x) <= bound for row in rows for x in row)
    lda, ldg = n + 3, m + 2
    flat = [value for row in rows for value in
            [*[math.ldexp(float(x), exponent) for x in row],
             SENTINEL, SENTINEL, SENTINEL]]
    source = (C.c_double * len(flat))(*flat)
    exact = [[sum(rows[i][k] * rows[j][k] for k in range(n))
              for j in range(i + 1)] for i in range(m)]
    results = []
    for function in functions:
        output = (C.c_double * (m * ldg))(*([SENTINEL] * (m * ldg)))
        status = function(source, m, n, lda, output, ldg, exponent, max_q)
        assert status == 0, (name, status)
        for i in range(m):
            for j in range(ldg):
                value = output[i * ldg + j]
                if j <= i:
                    assert Fraction.from_float(value) == \
                           exact_scale(exact[i][j], 2 * exponent), \
                           (name, i, j, value, exact[i][j])
                else:
                    assert math.isnan(value), (name, "upper/padding changed", i, j)
        results.append([output[i * ldg + j] for i in range(m) for j in range(i + 1)])
    assert results[0] == results[1], (name, "native/scalar disagreement")
    return {"name": name, "m": m, "n": n, "scale_exp": exponent,
            "max_q": max_q, "coefficient_cap": bound,
            "entries_per_implementation": m * (m + 1) // 2}


def fma_function():
    if hasattr(math, "fma"):
        return math.fma
    library = C.CDLL(C.util.find_library("m") or None)
    function = library.fma
    function.argtypes = [C.c_double, C.c_double, C.c_double]
    function.restype = C.c_double
    return function


def strict_counterexample():
    """Use actual correctly rounded FMA, then check each result as a Fraction."""
    fma = fma_function()
    b = 1 << 26
    d = 2 * b - 1
    start = fma(float(b), float(b), float(b * b))
    grouped_prefix = fma(float(d), float(b - 1), start)
    grouped_final = fma(float(-b), float(d), grouped_prefix)
    sequential_first_group = fma(float(d), float(-b), start)
    sequential_final = fma(float(b - 1), float(d), sequential_first_group)
    exact_prefix = 4 * b * b - 3 * b + 1
    exact_final = 2 * b * b - 2 * b + 1
    assert 2 * b * b == LIMIT
    assert Fraction.from_float(grouped_prefix) == exact_prefix - 1
    assert Fraction.from_float(grouped_final) == exact_final - 1
    assert Fraction.from_float(sequential_first_group) == b
    assert Fraction.from_float(sequential_final) == exact_final
    return {"rows": [[0, -b], [b - 1, b]], "n": 2, "q": 1, "B": b,
            "exact_grouped_prefix": exact_prefix,
            "rounded_grouped_prefix": int(grouped_prefix),
            "exact_diagonal": exact_final,
            "grouped_diagonal": int(grouped_final),
            "sequential_first_group": int(sequential_first_group),
            "sequential_diagonal": int(sequential_final)}


def finite_invariant_certificate():
    """Exhaust every n=2/B=1 state and all four support relations."""
    checked = 0
    support_relations = set()
    for hi, hj, ui, uj, ai, aj, ti, tj in itertools.product((-1, 0, 1), repeat=8):
        old, new, tail = ui * uj, ai * aj, ti * tj
        value = hi * hj + old
        assert abs(value) <= 2
        value -= hi * hj
        assert value == old and abs(value) <= 1
        value += tail
        assert value == old + tail and abs(value) <= 2
        di, dj = ai - ui, aj - uj
        support_relations.add((di != 0, dj != 0))
        operations = 2
        if di:
            value += di * uj
            operations += 1
        assert value == ai * uj + tail and abs(value) <= 2
        if dj:
            value += ai * dj
            operations += 1
        assert value == new + tail and abs(value) <= 2
        assert operations == 2 + int(di != 0) + int(dj != 0)
        checked += 1
    assert len(support_relations) == 4
    return {"exhaustive_n2_B1_states": checked,
            "support_relations": len(support_relations),
            "every_prefix_and_product_count": "passed",
            "scope": "finite supplemental certificate; arbitrary-n proof is analytical"}


def csr_group_certificate():
    """Check each CSR-group prefix against an independently updated factor list.

    An ordered-word polynomial check separately certifies the coordinate identity.
    The numeric cases include disjoint and overlapping supports at multiple n.
    """
    from collections import Counter
    def add(*terms):
        out = Counter()
        for term in terms:
            out.update(term)
        return {k: v for k, v in out.items() if v}
    def multiply(left, right):
        out = Counter()
        for a, x in left.items():
            for b, y in right.items():
                out[a + b] += x * y
        return dict(out)
    ui, uj = {("Ui",): 1}, {("Uj",): 1}
    di, dj = {("Di",): 1}, {("Dj",): 1}
    ai, aj = add(ui, di), add(uj, dj)
    lhs = add(multiply(ui, uj), multiply(di, uj), multiply(ai, dj))
    rhs = multiply(ai, aj)
    assert lhs == rhs and len(lhs) == 4
    rng = random.Random(2026101903)
    checked, prefixes = 0, 0
    relation_counts = {"unique_left": 0, "unique_right": 0, "overlap": 0}
    for n in [2, 3, 4, 7, 9, 17]:
        for pattern in range(16):
            old_i = [rng.randint(-3, 3) for _ in range(n - 1)]
            old_j = [rng.randint(-3, 3) for _ in range(n - 1)]
            new_i, new_j = old_i[:], old_j[:]
            for k in range(n - 1):
                if rng.randrange(4) < pattern % 4:
                    new_i[k] = rng.randint(-3, 3)
                if rng.randrange(4) < pattern // 4:
                    new_j[k] = rng.randint(-3, 3)
            left = [(k, new_i[k] - old_i[k]) for k in range(n - 1)
                    if new_i[k] != old_i[k]]
            right = [(k, new_j[k] - old_j[k]) for k in range(n - 1)
                     if new_j[k] != old_j[k]]
            si, sj = {k for k, _ in left}, {k for k, _ in right}
            relation_counts["unique_left"] += len(si - sj)
            relation_counts["unique_right"] += len(sj - si)
            relation_counts["overlap"] += len(si & sj)
            tail = rng.randint(-3, 3) * rng.randint(-3, 3)
            value = sum(a * b for a, b in zip(old_i, old_j)) + tail
            mixed_i = old_i[:]
            for k, d in left:
                value += d * old_j[k]
                mixed_i[k] = new_i[k]
                assert value == sum(a * b for a, b in zip(mixed_i, old_j)) + tail
                assert abs(value) <= n * 9
                prefixes += 1
            assert mixed_i == new_i
            mixed_j = old_j[:]
            for k, d in right:
                value += new_i[k] * d
                mixed_j[k] = new_j[k]
                assert value == sum(a * b for a, b in zip(new_i, mixed_j)) + tail
                assert abs(value) <= n * 9
                prefixes += 1
            assert mixed_j == new_j
            assert value == sum(a * b for a, b in zip(new_i, new_j)) + tail
            checked += 1
    assert all(relation_counts.values())
    return {"ordered_word_identity": "passed", "monomials": len(lhs),
            "csr_group_cases": checked, "csr_group_prefixes": prefixes,
            "support_coordinates": relation_counts,
            "scope": "finite symbolic/prefix checks; arbitrary-n theorem is analytical"}


def check_rejections_and_boundaries(functions):
    checks = set()
    one = (C.c_double * 1)(0.)
    for function in functions:
        output = (C.c_double * 36)(*([SENTINEL] * 36))
        assert function(None, 0, 2**63, 0, None, 0, -999, 999) == 0
        checks.add("m=0 unconditional no-op")
        assert function(None, 4, 0, 0, output, 6, 200, 64) == 0
        for i in range(4):
            for j in range(6):
                assert output[i * 6 + j] == 0 if j <= i else math.isnan(output[i * 6 + j])
        checks.add("n=0 lower zeros; upper/padding preserved")
        for n in [1, 2, 3, 7, 8, 9, 17, 65, 65536]:
            source = (C.c_double * n)(*([0.] * n))
            source[0] = float(cap(n) + 1)
            assert function(source, 1, n, n, output, 1, 0, 0) == -4, n
            source[0] = -source[0]
            assert function(source, 1, n, n, output, 1, 0, 0) == -4, n
        checks.add("positive/negative first coefficient beyond integer-adjusted cap")
        for value, exponent, code, label in [
            (math.inf, 0, -2, "+infinity"), (-math.inf, 0, -2, "-infinity"),
            (math.nan, 0, -2, "NaN"), (0.5, 0, -4, "half lattice step"),
            (math.ldexp(1., -1000), 200, -4, "tiny input scales to zero"),
            (math.ulp(0.), 200, -4, "smallest subnormal"),
            (math.ldexp(1., -401), -400, -4, "half step at exponent minimum"),
            (math.nextafter(math.ldexp(1., 200), math.inf), 200, -4,
             "one ulp off grid"),
        ]:
            for n in [1, 8, 9, 16, 17, 24, 25]:
                for lane in range(n):
                    source = (C.c_double * n)(*([0.] * n))
                    source[lane] = value
                    assert function(source, 1, n, n, output, 1, exponent, 0) == code, \
                           (label, n, lane)
            checks.add(label)
        for exponent in [-401, 201]:
            assert function(one, 1, 1, 1, output, 1, exponent, 0) == -4
        assert function(one, 1, 65537, 65537, output, 1, 0, 0) == -4
        assert function(one, 1, 1, 1, output, 1, 0, 65) == -4
        assert function(None, 1, 1, 1, output, 1, 0, 0) == -1
        assert function(one, 1, 1, 0, output, 1, 0, 0) == -1
        assert function(one, 1, 1, 1, None, 1, 0, 0) == -1
        assert function(one, 2, 1, 1, output, 1, 0, 0) == -1
        size_max = C.c_size_t(-1).value
        assert function(one, size_max, 1, 1, output, size_max, 0, 0) == -1
        checks.add("invalid pointers, strides, extents, exponent/n/q limits")
        source = (C.c_double * 4)(1., 2., 3., 4.)
        assert function(source, 2, 2, 2, output, 2, 0, 0) == -2
        checks.add("actual defects exceed caller cap")
        zeros = (C.c_double * 8)(*([-0., +0.] * 4))
        assert function(zeros, 1, 8, 8, output, 1, 200, 0) == 0 and output[0] == 0
        checks.add("both zero signs admitted")
    return sorted(checks)


def main():
    assert hashlib.sha256(ORIGINAL.read_bytes()).hexdigest() == ORIGINAL_SHA256, \
           "The copied scaffold changed; provenance must be reviewed before comparison."
    assert hashlib.sha256(SUBSTITUTE.read_bytes()).hexdigest() == SUBSTITUTE_SHA256, \
           "The coordinate-substitution source must remain frozen."
    builds = [build(False), build(True)]
    functions = [item[0] for item in builds]
    rng = random.Random(2026101902)
    cases = []
    for n in range(1, 26):
        q = min(3, n - 1)
        cases.append((f"all_vector_tails_n{n}", shifted(rng, 7, n, q, cap(n)),
                      [-400, -1, 0, 200][n % 4], q))
    for n in [1, 2, 3, 7, 8, 9, 17, 65, 65536]:
        cases.append((f"maximal_uniform_anchor_n{n}", [[cap(n)] * n], -400, 64))
    for exponent in [-400, -399, -1, 0, 1, 199, 200]:
        cases.append((f"q0_exponent_{exponent}", shifted(rng, 9, 11, 0, cap(11)),
                      exponent, 0))
    b = 1 << 26
    cases += [
        ("strict_grouped_failure", [[0, -b], [b - 1, b]], 0, 1),
        ("maximum_residual_2B", [[-cap(9)] * 9, [cap(9)] * 9], 200, 8),
        ("maximum_q64", shifted(rng, 4, 65, 64, cap(65)), -400, 64),
        ("recurrence_chain_m37", shifted(rng, 37, 3, 1, cap(3)), 0, 1),
        ("zero_rows", [[0] * 9 for _ in range(7)], -400, 0),
        ("n1_signed_extremes", [[-cap(1)], [cap(1)], [0], [1], [-1]], 200, 0),
    ]
    # Exercise each shifted lane separately, especially both ends of each
    # full 8-value block and every reserved scalar tail.
    for n in [9, 16, 17, 24, 25]:
        for k in range(n - 1):
            a0 = list(range(n))
            a1 = a0[1:] + [n]
            a1[k] = -a1[k]
            cases.append((f"single_defect_n{n}_k{k}", [a0, a1], 200, 1))
    reports = [check_case(functions, *case) for case in cases]
    report = {
        "status": "all sequential exact comparisons and rejection checks passed",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "original_source_sha256": ORIGINAL_SHA256,
        "scaffold_source_sha256": SUBSTITUTE_SHA256,
        "timing": "none: correctness only; no timers or vendor calls",
        "reference": "direct Python integer Gram and exact Fraction scaling",
        "implementations": ["native (NEON on arm64)", "forced scalar"],
        "compiler_commands": [item[1] for item in builds],
        "strict_counterexample": strict_counterexample(),
        "finite_invariant_certificate": finite_invariant_certificate(),
        "csr_group_certificate": csr_group_certificate(),
        "cases": reports,
        "entries_per_implementation": sum(x["entries_per_implementation"] for x in reports),
        "rejection_and_boundary_classes": check_rejections_and_boundaries(functions),
    }
    (BASE / "verification/sequential_results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ["status", "source_sha256",
                     "entries_per_implementation", "timing"]}))


if __name__ == "__main__":
    main()
