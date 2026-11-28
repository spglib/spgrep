# Algorithm for enumerating irreps

Irreps of crystallographic point groups or space groups are calculated in following steps:

Method-A ([`spgrep.irreps.enumerate_unitary_irreps_from_regular_representation`][spgrep.rep.enumerate.enumerate_unitary_irreps_from_regular_representation])

1. Construct regular representation (crystallographic point group) or projective regular representation (space group).

2. Construct random matrix that commutes with (projective) regular representation and obtain irreps by diagonalizing it.

Method-B ([`spgrep.irreps.enumerate_unitary_irreps_from_solvable_group_chain`][spgrep.symmetry.enumerate.enumerate_unitary_irreps_from_solvable_group_chain])

1. Decompose little co-group $\overline{\mathcal{G}}^{\mathbf{k}}$ into series

    $$
    \begin{aligned}
    1 = G_{0} \triangleleft G_{1} \triangleleft \dots \triangleleft G_{m} = \overline{\mathcal{G}}^{\mathbf{k}}.
    \end{aligned}
    $$

2. Construct irrep $\Delta$ of $G_{i}$ and use induced representation $\Delta \uparrow G_{i+1}$ to obtain irreps of $G_{i+1}$

## Subpages

- [On-the-fly irreps generation from regular representation](irreps_from_regular.md)
- [On-the-fly irreps generation from solvable group chain](irreps_from_chain.md)
- [Intertwiner](intertwiner.md)
