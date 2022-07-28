# spgrep
[![testing](https://github.com/spglib/spgrep/actions/workflows/testing.yml/badge.svg)](https://github.com/spglib/spgrep/actions/workflows/testing.yml)

**spgrep** is a Python package of on-the-fly generator of space-group irreducible representations.

## Features

- Calculate irreducible representations (irreps) of space groups from spglib’s cell and kpoints
- Calculate irreps of crystallographic point groups
- Find symmetry-adapted basis forming given irreps
- Minimal dependencies (numpy and spglib)

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

Installation
```shell
conda create -n spgrep python=3.10 pip
conda activate spgrep
git clone git@github.com:spglib/spgrep.git
cd spgrep
pip install -e ".[dev,docs]"
pre-commit install
```

Document
```shell
sphinx-autobuild docs docs_build
# open localhost:8000 in your browser
```

Confirm the version number via `setuptools-scm`
```shell
python -m setuptools_scm
# e.g. 0.1.dev1+g7aa164a.d20220704
```
