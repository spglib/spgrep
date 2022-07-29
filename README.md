# spgrep
[![testing](https://github.com/spglib/spgrep/actions/workflows/testing.yml/badge.svg)](https://github.com/spglib/spgrep/actions/workflows/testing.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/spglib/spgrep/main.svg)](https://results.pre-commit.ci/latest/github/spglib/spgrep/main)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spgrep)
[![PyPI version](https://badge.fury.io/py/spgrep.svg)](https://badge.fury.io/py/spgrep)

**spgrep** is a Python package of on-the-fly generator of space-group irreducible representations.

- Github: https://github.com/spglib/spgrep
- Document: https://spglib.github.io/spgrep/
- PyPI: https://pypi.org/project/spgrep

## Features

- Calculate irreducible representations (irreps) of space groups from spglib’s cell and kpoints
- Calculate irreps of crystallographic point groups
- Calculate physically irreducible representations (irreps over real numbers)
- Find symmetry-adapted basis forming given irreps
- Minimal dependencies (numpy and spglib)

## Usage

```python
from spgrep import get_spacegroup_irreps
from spgrep.representation import get_character

# Rutile structure (https://materialsproject.org/materials/mp-2657/)
# P4_2/mnm (No. 136)
a = 4.603
c = 2.969
x_4f = 0.3046
lattice = [
    [a, 0, 0],
    [0, a, 0],
    [0, 0, c],
]
positions = [
    [0, 0, 0],  # Ti(2a)
    [0.5, 0.5, 0.5],  # Ti(2a)
    [x_4f, x_4f, 0],  # O(4f)
    [-x_4f, -x_4f, 0],  # O(4f)
    [-x_4f + 0.5, x_4f + 0.5, 0.5],  # O(4f)
    [x_4f + 0.5, -x_4f + 0.5, 0.5],  # O(4f)
]
numbers = [0, 0, 1, 1, 1, 1]

kpoint = [0.5, 0, 0]  # X point
irreps, rotations, translations, mapping_little_group = get_spacegroup_irreps(
    lattice, positions, numbers, kpoint
)

# Symmetry operations by spglib
assert len(rotations) == 16
assert len(translations) == 16

# At X point, the little co-group is isomorphic to mmm (order=8)
assert len(mapping_little_group) == 8
print(mapping_little_group)  # [ 0,  1,  4,  5,  8,  9, 12, 13]

# Two two-dimensional irreps
for irrep in irreps:
    print(get_character(irrep))
# [2.+0.j 0.+0.j 0.+0.j 2.+0.j 0.+0.j 0.+0.j 0.+0.j 0.+0.j]
# [2.+0.j 0.+0.j 0.+0.j -2.+0.j 0.+0.j 0.+0.j 0.+0.j 0.+0.j]
```

## Installation

```shell
conda create -n spgrep python=3.10 pip
conda activate spgrep
git clone git@github.com:spglib/spgrep.git
cd spgrep
pip install -e .
```

## License

spgrep is released under a BSD 3-clause license.

## Development

### Installation

```shell
conda create -n spgrep python=3.10 pip
conda activate spgrep
git clone git@github.com:spglib/spgrep.git
cd spgrep
pip install -e ".[dev,docs]"
pre-commit install

# Run pre-commit manually
pre-commit run --all-file 
```

### Document

```shell
sphinx-autobuild docs docs_build
# open localhost:8000 in your browser
```

### Release

```shell
# Confirm the version number via `setuptools-scm`
python -m setuptools_scm
# Test a release version at rc branch
git checkout rc

# Update changelog here

# Push with tag
git tag <next-version>
git push --follow-tags origin rc

# Merge to main
git checkout main
git merge rc
git push --follow-tags origin main
```
