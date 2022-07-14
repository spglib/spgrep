# Formulation

General references [^BC09] [^ITO96]

## Overview

Irreps of crystallographic point groups or space groups are calculated in following steps:

Option-A
1. Construct regular representation (crystallographic point group) or projective regular representation (space group).
2. Construct random matrix that commutes with (projective) regular representation and obtain irreps by diagonalizing it.

Option-B
1. Decompose little co-group {math}`\overline{\mathcal{G}}^{\mathbf{k}}` into series
  ```{math}
    E = G_{0} \triangleleft G_{1} \triangleleft \dots \triangleleft G_{m} = \overline{\mathcal{G}}^{\mathbf{k}}.
  ```
2. Use induced representation


```{toctree}
  Irreps of space group <spacegroup_irreps>
  On-the-fly irreps generation from regular representation <irreps_from_regular>
  On-the-fly irreps generation from solvable group chain <irreps_generation>
```

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[^BC09]: C. Bradley and A. P. Cracknell, The mathematical theory of symmetry in solids (Oxford, London, 2009).
[^ITO96]: T. Inui, Y. Tanabe, and Y. Onodera, Group theory and its applications in physics (Springer, Berlin, 1996).
