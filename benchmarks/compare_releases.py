"""Compare two committed revisions in one Python environment.

    python benchmarks/compare_releases.py --baseline v0.7.0 --candidate HEAD \
        --output benchmarks/results/v0.8.0.json

Uses the current benchmark_cayley_table.py for both revisions. Source archives
are temporary directories, not Git worktrees. Each batch uses a fresh process;
imports, input construction, validation, and warmup are outside the timed region.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from statistics import median


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default="v0.7.0")
    parser.add_argument("--candidate", default="HEAD")
    parser.add_argument("--repeat", type=int, default=9)
    parser.add_argument("--number", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.repeat < 1 or args.number < 1:
        parser.error("--repeat and --number must be positive")

    root = Path(__file__).resolve().parents[1]
    benchmark = Path(__file__).with_name("benchmark_cayley_table.py").resolve()
    revisions = {
        name: subprocess.check_output(
            ["git", "rev-parse", f"{ref}^{{commit}}"], cwd=root, text=True
        ).strip()
        for name, ref in (("baseline", args.baseline), ("candidate", args.candidate))
    }
    environment = os.environ.copy()
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[name] = "1"

    batches: dict[str, list[dict]] = {name: [] for name in revisions}
    with tempfile.TemporaryDirectory(prefix="spgrep-release-") as directory:
        snapshots = {}
        for name, revision in revisions.items():
            snapshot = Path(directory) / name
            archive = subprocess.check_output(
                ["git", "archive", "--format=zip", revision, "src"], cwd=root
            )
            with zipfile.ZipFile(io.BytesIO(archive)) as source:
                source.extractall(snapshot)
            snapshots[name] = snapshot

        def run(name, number):
            output = Path(directory) / f"{name}.json"
            environment["PYTHONPATH"] = str(snapshots[name] / "src")
            subprocess.run(
                [
                    sys.executable,
                    str(benchmark),
                    "--repeat",
                    "1",
                    "--number",
                    str(number),
                    "--output",
                    str(output),
                ],
                cwd=root,
                env=environment,
                check=True,
                stdout=subprocess.DEVNULL,
            )
            result = json.loads(output.read_text())
            imported = Path(result["metadata"]["source"]).resolve()
            assert imported == (snapshots[name] / "src/spgrep/__init__.py").resolve()
            result["metadata"]["source"] = "src/spgrep/__init__.py"
            return result

        probes = {name: run(name, 1) for name in revisions}
        estimate = (
            sum(case["median_ms"] for probe in probes.values() for case in probe["results"])
            * args.number
            * args.repeat
            / 1000
        )
        print(f"Estimated timed work: {estimate:.1f} s, plus process setup and warmup", flush=True)
        for index in range(args.repeat):
            order = list(revisions)
            if index % 2:
                order.reverse()
            for name in order:
                batches[name].append(run(name, args.number))
            print(f"Completed paired batch {index + 1}/{args.repeat}", flush=True)

    metadata = batches["baseline"][0]["metadata"]
    assert all(batch["metadata"] == metadata for values in batches.values() for batch in values)
    cases = [case["case"] for case in batches["baseline"][0]["results"]]
    assert all(
        [case["case"] for case in batch["results"]] == cases
        for values in batches.values()
        for batch in values
    )
    results = []
    for case in cases:
        samples = {
            name: [
                next(row for row in batch["results"] if row["case"] == case)["median_ms"]
                for batch in values
            ]
            for name, values in batches.items()
        }
        before, after = (median(samples[name]) for name in revisions)
        results.append(
            {
                "case": case,
                "baseline_ms": before,
                "candidate_ms": after,
                "speedup": before / after,
                "samples_ms": samples,
            }
        )
        print(f"{case:<22} {before:9.3f} -> {after:9.3f} ms  {before / after:6.2f}x")
    metadata.update(
        {
            "measured_at": datetime.now(timezone.utc).isoformat(),
            "revisions": revisions,
            "repeat": args.repeat,
            "benchmark_sha256": hashlib.sha256(benchmark.read_bytes()).hexdigest(),
            "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"metadata": metadata, "results": results}, indent=2) + "\n")


if __name__ == "__main__":
    main()
