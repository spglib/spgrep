# Formulation

General references {cite}`Bradley2009-ze,Inui1996-et`

## Overview

Irreps of crystallographic point groups or space groups are calculated in following steps:

Option-A ({func}`spgrep.irreps.enumerate_unitary_irreps_from_regular_representation`)
1. Construct regular representation (crystallographic point group) or projective regular representation (space group).
2. Construct random matrix that commutes with (projective) regular representation and obtain irreps by diagonalizing it.

Option-B ({func}`spgrep.irreps.enumerate_unitary_irreps_from_solvable_group_chain`)
1. Decompose little co-group {math}`\overline{\mathcal{G}}^{\mathbf{k}}` into series
  ```{math}
    E = G_{0} \triangleleft G_{1} \triangleleft \dots \triangleleft G_{m} = \overline{\mathcal{G}}^{\mathbf{k}}.
  ```
2. Construct irrep {math}`\Delta` of {math}`G_{i}` and use induced representation {math}`\Delta \uparrow G_{i+1}` to obtain irreps of {math}`G_{i+1}`

## Subpages

```{toctree}
  :maxdepth: 1
  Irreps of space group <spacegroup_irreps>
  On-the-fly irreps generation from regular representation <irreps_from_regular>
  On-the-fly irreps generation from solvable group chain <irreps_from_chain>
  Reality of irreps <reality>
  Symmetry-adapted basis <projection>
  Intertwiner <intertwiner>
  Spin representation <spinor>
```

## References

```{bibliography}
:filter: docname in docnames
```
