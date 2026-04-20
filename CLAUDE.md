# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`spgrep` is a Python library that generates irreducible representations (irreps) of crystallographic point groups and space groups on the fly from spglib symmetry data. It supports linear irreps, physically (real) irreps, projective spinor irreps, and projective irreducible co-representations (with antiunitary / time-reversal operations).

Python 3.10+. Runtime deps are intentionally minimal: `numpy`, `spglib`, `hsnf`, `typing-extensions`.

## Commands

The project uses `uv` for environment management and `prek` (a pre-commit reimpl) for hooks. A `justfile` wraps the common ones.

```shell
just install         # uv sync --all-extras
just test            # uv run pytest -v
just prek            # prek run --all-files (ruff-format, ruff, mypy, uv-lock, ...)
just docs            # uv run sphinx-autobuild docs docs_build
```

Direct equivalents (use these if `just` is unavailable):

```shell
uv sync --all-extras
uv run pytest -v
uv run pytest tests/test_core.py::test_get_spacegroup_irreps -v   # single test
uv run pytest -v --cov=spgrep --cov-report=xml tests/             # CI invocation
prek run --all-files                                              # lint + format + mypy
```

CI (`.github/workflows/testing.yml`) runs the matrix over Python 3.10/3.11/3.12/3.13 and a separate docs build job; pandoc is required for docs builds.

Releases are tag-driven: bump `docs/changelog.md`, then `git tag <version> && git push origin <version>`. Version is derived by `setuptools-scm`.

## Architecture

`spgrep` is layered. Lower layers are crystal-agnostic group/representation code; upper layers add space-group/little-group/spin specifics; `core.py` is the user-facing facade.

```
spgrep.core              <- public API (re-exported from spgrep/__init__.py)
   |
   +-- spgrep.corep              co-representations (antiunitary ops, magnetic groups)
   +-- spgrep.spinor             spinor factor systems and spin small reps
   +-- spgrep.symmetry/*         space-group / little-group machinery
          |
          +-- spgrep.rep/*       generic finite-group representation toolkit
                                 (no crystal knowledge)
```

### `spgrep.rep/` — generic finite-group toolkit
Pure abstract group/representation primitives. No notion of crystals, k-points, or rotations.
- `group.py` — Cayley table indexing, identity/inverse lookup, group order
- `representation.py` — characters, intertwiners, direct products, representation checks
- `irreps.py` — Frobenius-Schur indicator, equivalence of irreps, projection to an irrep
- `enumerate.py` — decompose a (projective) regular representation into unitary irreps; this is the algorithmic core (`enumerate_unitary_irreps_from_regular_representation`)
- `pir.py` — physically (real) irreducible representations

### `spgrep.symmetry/` — space-group layer
Builds on `rep/` by feeding in space-group data.
- `group.py` — little group construction, factor systems from little group
- `representation.py` — projective regular representation derived from a factor system
- `enumerate.py` — `enumerate_small_representations`, `enumerate_unitary_irreps`, irrep purification
- `pointgroup.py` — point-group chain generators (used by the "Neto" method)
- `transform.py` — primitive-cell transformations and unique primitive symmetry filtering
- `subgroup.py`, `tensor.py` — additional space-group utilities

### `spgrep.core`, `spgrep.corep`, `spgrep.spinor` — user-facing facade
- `core.py` exports the six top-level entry points (`get_spacegroup_irreps`, `get_spacegroup_irreps_from_primitive_symmetry`, the `_spinor_` variants, and the `_crystallographic_pointgroup_` variants). These are re-exported from `spgrep/__init__.py`.
- `corep.py` handles co-representations (antiunitary operations).
- `spinor.py` provides spin factor systems and small spin representations.

### Two algorithm choices
Most enumeration entry points take `method: Literal["Neto", "random"]`:
- `"Neto"` — construct irreps along a fixed chain of subgroups of the little co-group (default, deterministic).
- `"random"` — diagonalize a random matrix that commutes with the regular representation. Controlled by `max_num_random_generations`.

### Deprecated top-level shims
`spgrep/representation.py`, `spgrep/group.py`, `spgrep/pointgroup.py`, `spgrep/transform.py`, `spgrep/tensor.py`, and `spgrep/irreps.py` are kept as deprecated re-export shims that emit `DeprecationWarning`. New code should import from `spgrep.rep.*` and `spgrep.symmetry.*` instead. When changing implementations, edit the new locations and let the shim forward.

### `spgrep._constants`
Default tolerances (`RTOL`, `ATOL`) and `MAX_NUM_RANDOM_GENERATIONS` live here and are threaded through every public API as keyword arguments.

## Tests

Tests in `tests/` mirror the source layout (`tests/rep/`, `tests/symmetry/`, top-level). `pyproject.toml` pins `testpaths = ["tests"]`. Fixtures are in `tests/conftest.py`.

## Conventions specific to this repo

- Tolerances: numerical results are compared with `RTOL`/`ATOL` from `spgrep._constants`. When adding numerical comparisons in new code or tests, prefer these constants over ad-hoc literals.
- All public APIs use `from __future__ import annotations` and NumPy-style docstrings; the array-shape conventions in those docstrings are load-bearing (they document index semantics callers rely on).
- The `_constants.MAX_NUM_RANDOM_GENERATIONS` retry knob exists because the `"random"` method can statistically fail; honor it rather than introducing infinite loops.
- Ruff is configured with `line-length = 99` and selects `E`, `F`, `I`, `UP`. `D2xx` pydocstyle rules listed in `[tool.ruff.lint.extend-ignore]` are intentionally disabled — do not re-enable them piecemeal.
- mypy excludes `docs/` and sets `warn_no_return = false`.
