# Formulation

## Overview

Irreps of crystallographic point groups or space groups are calculated in following steps:

1. Construct regular representation (crystallographic point group) or projective regular representation (space group).
1. Construct random matrix that commutes with (projective) regular representation and obtain irreps by diagonalizing it.
1. Symmetrize the obtained irrep matrices by subduced representation.

```{toctree}
  Little group <little_group>
  On-the-fly irreps generation <irreps_generation>
```
