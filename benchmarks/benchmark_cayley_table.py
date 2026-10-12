"""Measure Cayley-table construction and public irrep calls.

Run before and after a change in the same environment:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
        PYTHONPATH=src python benchmarks/benchmark_cayley_table.py \
        --output .cache/cayley.json

Input construction, imports, and one warmup call are outside the timed region.
Public-call timings include symmetry search or spinor preparation as applicable.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
from functools import partial
from pathlib import Path
from statistics import median
from timeit import repeat

import numpy as np

from spgrep import get_spacegroup_irreps
from spgrep.core import get_crystallographic_pointgroup_spinor_irreps_from_symmetry
from spgrep.symmetry.group import get_cayley_table
from spgrep.symmetry.pointgroup import pg_dataset

CASES = ("table-C3v", "table-Oh", "table-gray-Oh", "crystal-Neto", "crystal-random", "spinor-Oh")


def make_case(name):
    rotations = np.array(pg_dataset["3m" if name == "table-C3v" else "m-3m"][0])
    if name == "table-gray-Oh":
        time_reversals = np.repeat([0, 1], len(rotations))
        rotations = np.concatenate([rotations, rotations])
        return partial(get_cayley_table, rotations, time_reversals), len(rotations)
    if name.startswith("table-"):
        return partial(get_cayley_table, rotations), len(rotations)
    if name == "spinor-Oh":
        return (
            partial(
                get_crystallographic_pointgroup_spinor_irreps_from_symmetry, np.eye(3), rotations
            ),
            len(rotations),
        )
    return (
        partial(
            get_spacegroup_irreps,
            np.eye(3),
            np.zeros((1, 3)),
            np.array([1]),
            np.zeros(3),
            method="random" if name == "crystal-random" else "Neto",
        ),
        len(rotations),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES, action="append")
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--number", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeat < 1 or args.number < 1:
        parser.error("--repeat and --number must be positive")

    metadata = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
        "thread_limits": {
            name: os.environ.get(name)
            for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
        },
        "repeat": args.repeat,
        "number": args.number,
    }
    print(json.dumps(metadata), flush=True)
    print(f"{'case':<18} {'order':>5} {'median ms':>12}", flush=True)
    results = []
    for name in args.case or CASES:
        call, order = make_case(name)
        result = call()
        if name.startswith("table-"):
            assert result.shape == (order, order)
            np.testing.assert_array_equal(
                np.sort(result, axis=1), np.broadcast_to(np.arange(order), (order, order))
            )
        else:
            irreps = result[0]
            assert sum(irrep.shape[1] ** 2 for irrep in irreps) == order
            characters = np.array([np.trace(irrep, axis1=1, axis2=2) for irrep in irreps])
            np.testing.assert_allclose(
                characters.conj() @ characters.T / order, np.eye(len(irreps)), atol=1e-8
            )
        samples_ms = [
            1000 * elapsed / args.number
            for elapsed in repeat(call, number=args.number, repeat=args.repeat)
        ]
        median_ms = median(samples_ms)
        results.append(
            {"case": name, "order": order, "median_ms": median_ms, "samples_ms": samples_ms}
        )
        print(f"{name:<18} {order:5d} {median_ms:12.3f}", flush=True)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps({"metadata": metadata, "results": results}, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
