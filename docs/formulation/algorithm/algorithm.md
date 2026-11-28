# Algorithm for enumerating irreps

Understand the two algorithms spgrep uses to enumerate irreps of crystallographic point groups and space groups.

!!! info "Before you start"

    Read [space-group irreps](../irreps/spacegroup_irreps.md), especially little groups
    and projective representations. For a runnable calculation, use the
    [quickstart](../../quickstart.md).

The core functions use the subgroup-chain method by default (`method="Neto"`).
The regular-representation method is selected with `method="random"`.
Both compute irreps up to equivalence; their basis choices and ordering can differ.
See the [core API reference](../../api/api_core.md) for numerical parameters.

## Regular-representation method

[`spgrep.rep.enumerate.enumerate_unitary_irreps_from_regular_representation`][spgrep.rep.enumerate.enumerate_unitary_irreps_from_regular_representation]

1. Construct regular representation (crystallographic point group) or projective regular representation (space group).

2. Construct random matrix that commutes with (projective) regular representation and obtain irreps by diagonalizing it.

## Subgroup-chain method

[`spgrep.symmetry.enumerate.enumerate_unitary_irreps_from_solvable_group_chain`][spgrep.symmetry.enumerate.enumerate_unitary_irreps_from_solvable_group_chain]

1. Decompose little co-group $\overline{\mathcal{G}}^{\mathbf{k}}$ into series

    $$
    \begin{aligned}
    1 = G_{0} \triangleleft G_{1} \triangleleft \dots \triangleleft G_{m} = \overline{\mathcal{G}}^{\mathbf{k}}.
    \end{aligned}
    $$

2. Construct irrep $\Delta$ of $G_{i}$ and use induced representation $\Delta \uparrow G_{i+1}$ to obtain irreps of $G_{i+1}$

## Next steps

- [On-the-fly irreps generation from regular representation](irreps_from_regular.md)
- [On-the-fly irreps generation from solvable group chain](irreps_from_chain.md)
- [Intertwiner](intertwiner.md)
