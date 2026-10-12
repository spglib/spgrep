"""Benchmark irrep extraction and its caller with deterministic inputs.

Run from the repository root before and after a change, using the same environment:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
        PYTHONPATH=src python benchmarks/benchmark_irreps_from_matrix.py

Timings are medians of repeated batches, excluding imports and input construction.
The caller timing includes random-matrix construction and decomposition. The helper
timing uses the exact same representation and commuting matrix, captured from the
caller. C64 stresses many one-dimensional irreps; dense-C3v exercises a dense,
non-regular representation whose dimension differs from the group order.
"""

from __future__ import annotations

import argparse
import platform
from functools import partial
from itertools import permutations
from statistics import median
from timeit import repeat
from unittest.mock import patch

import numpy as np

import spgrep.rep.enumerate as enumeration
from spgrep.spinor import get_spinor_factor_system
from spgrep.symmetry.pointgroup import pg_dataset
from spgrep.symmetry.representation import get_projective_regular_representation

CASES = ("C3v", "Oh", "spinor-Oh", "C64", "dense-C3v")


def make_case(name):
    if name == "dense-C3v":
        rep = np.zeros((6, 4, 4), dtype=complex)
        for i, permutation in enumerate(permutations(range(3))):
            rep[i, :3, :3] = np.eye(3)[:, permutation]
            rep[i, 3, 3] = np.linalg.det(rep[i, :3, :3])
        rng = np.random.default_rng(0)
        basis, _ = np.linalg.qr(rng.random((4, 4)) + 1j * rng.random((4, 4)))
        rep = basis.conj().T @ rep @ basis
        return rep, enumeration.decompose_representation, [1, 1, 2]

    if name == "C64":
        order = 64
        rep = np.zeros((order, order, order), dtype=complex)
        for i in range(order):
            rep[i, (i + np.arange(order)) % order, np.arange(order)] = 1
        dimensions = [1] * order
    else:
        rotations = np.array(pg_dataset["3m" if name == "C3v" else "m-3m"][0])
        factor_system = np.ones((len(rotations), len(rotations)), dtype=complex)
        if name == "spinor-Oh":
            factor_system, _ = get_spinor_factor_system(np.eye(3), rotations)
        rep = get_projective_regular_representation(rotations, factor_system)
        dimensions = {
            "C3v": [1, 1, 2],
            "Oh": [1, 1, 1, 1, 2, 2, 3, 3, 3, 3],
            "spinor-Oh": [2, 2, 2, 2, 4, 4],
        }[name]
    return rep, enumeration.enumerate_unitary_irreps_from_regular_representation, dimensions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES, action="append")
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--number", type=int, default=20)
    args = parser.parse_args()
    if args.repeat < 1 or args.number < 1:
        parser.error("--repeat and --number must be positive")

    print(f"Python {platform.python_version()}, NumPy {np.__version__}, {platform.platform()}")
    print(f"Median of {args.repeat} batches of {args.number} calls; milliseconds per call")
    print(f"{'case':<12} {'order':>5} {'dim':>5} {'helper':>12} {'caller':>12}")
    for name in args.case or CASES:
        rep, caller, dimensions = make_case(name)
        with patch.object(
            enumeration, "_get_irreps_from_matrix", wraps=enumeration._get_irreps_from_matrix
        ) as capture:
            irreps = caller(rep)
        assert [irrep.shape[1] for irrep in irreps] == dimensions
        characters = np.array([np.trace(irrep, axis1=1, axis2=2) for irrep in irreps])
        np.testing.assert_allclose(
            characters.conj() @ characters.T / len(rep), np.eye(len(irreps)), atol=1e-8
        )
        assert capture.call_count == 1, "This case needs retries; choose a representative matrix"
        helper_args, helper_kwargs = capture.call_args
        helper = partial(enumeration._get_irreps_from_matrix, *helper_args, **helper_kwargs)
        helper_ms = (
            1000 * median(repeat(helper, number=args.number, repeat=args.repeat)) / args.number
        )
        caller_ms = (
            1000
            * median(repeat(partial(caller, rep), number=args.number, repeat=args.repeat))
            / args.number
        )
        print(
            f"{name:<12} {rep.shape[0]:5d} {rep.shape[1]:5d} {helper_ms:12.3f} {caller_ms:12.3f}"
        )


if __name__ == "__main__":
    main()
