# Spgrep

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/spglib/spgrep/develop/docs/_static/logo-dark.svg">
    <img alt="spgrep" src="https://raw.githubusercontent.com/spglib/spgrep/develop/docs/_static/logo.svg" width="360">
  </picture>
</p>

[![testing](https://github.com/spglib/spgrep/actions/workflows/testing.yml/badge.svg)](https://github.com/spglib/spgrep/actions/workflows/testing.yml)
[![codecov](https://codecov.io/gh/spglib/spgrep/graph/badge.svg?token=DQGVFCTB1P)](https://codecov.io/gh/spglib/spgrep)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spgrep)](https://img.shields.io/pypi/pyversions/spgrep)
[![PyPI version](https://badge.fury.io/py/spgrep.svg)](https://badge.fury.io/py/spgrep)
[![PyPI Downloads](https://img.shields.io/pypi/dm/spgrep)](https://img.shields.io/pypi/dm/spgrep)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.05269/status.svg)](https://doi.org/10.21105/joss.05269)

spgrep is a Python package for calculating irreducible representations of space groups
and crystallographic point groups from crystal structures or symmetry operations.

> **Who this is for:** Researchers who use crystal symmetry and Python arrays.
> Start with the [quickstart](docs/quickstart.md) if you are new to representation calculations.
> The [theory pages](docs/formulation/formulation.md) explain the mathematics separately.

<a id="installation"></a>
<a id="usage"></a>

## Start a calculation

Install the package with Python 3.10 or later:

```shell
python -m pip install spgrep
```

Follow the [quickstart](docs/quickstart.md) to calculate the irreducible representations
of rutile at the X point and associate each matrix with a symmetry operation.
Then read [how spgrep represents symmetry](docs/conventions.md) before using your own data.
For a source checkout and contributor tools, see [development](docs/development/development.md).

## Find the guide for your task

- To calculate irreps from operations you already have, use
  [Calculate irreps from primitive symmetry](docs/howto/primitive-symmetry.md).
- To construct tensor bases, analyze lattice vibrations, or work with spinors, choose
  a [worked tutorial](docs/examples/examples.md). You can skip the derivations to run these examples.
- To understand little groups, real representations, or spinor factor systems, read
  the [theory guide](docs/formulation/formulation.md). Start with space-group irreps,
  then choose the extension your calculation needs.
- To look up arguments, return arrays, or numerical options, use the
  [API reference](docs/api/api.md). Experienced users can start there.

## Supported calculations

spgrep enumerates complex irreducible representations, physically irreducible
representations over real vector spaces, spinor representations, and spinor
co-representations for magnetic symmetry. It also constructs symmetry-adapted bases.
See [core functions](docs/api/api_core.md) for the supported combinations of input
and representation type.

The [source repository](https://github.com/spglib/spgrep),
[documentation site](https://spglib.github.io/spgrep/), and
[PyPI package](https://pypi.org/project/spgrep/) provide the code, guides, and releases.

## Change log

See the [change log](docs/changelog.md) for recent changes.

## How to cite spgrep

If you use spgrep in your research, please consider citing the following work:

- [Kohei Shinohara, Atsushi Togo, Isao Tanaka, spgrep: On-the-fly generator of space-group irreducible representations, J. Open Source Softw. 8, 85, 5269 (2023)](https://doi.org/10.21105/joss.05269)

```
@article{spgrep,
    doi = {10.21105/joss.05269},
    url = {https://doi.org/10.21105/joss.05269},
    year = {2023},
    publisher = {The Open Journal},
    volume = {8},
    number = {85},
    pages = {5269},
    author = {Shinohara, Kohei and Togo, Atsushi and Tanaka, Isao},
    title = {spgrep: On-the-fly generator of space-group irreducible representations},
    journal = {J. Open Source Softw.}
}
```

## How to contribute

We welcome any contributions to improve functionalities.
Please open [issues](https://github.com/spglib/spgrep/issues?q=is%3Aissue+is%3Aopen+sort%3Aupdated-desc) or create [pull requests](https://github.com/spglib/spgrep/pulls?q=is%3Apr+is%3Aopen+sort%3Aupdated-desc).
See the [development page](docs/development/development.md) for preparing an environment.

## License

Spgrep is released under a BSD 3-clause license.
