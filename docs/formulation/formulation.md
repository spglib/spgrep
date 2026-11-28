# Understand the representations and algorithms

Choose the mathematical explanation that supports your representation calculation.

!!! info "Before you start"

    These pages assume linear algebra and crystallographic symmetry.
    Read [how spgrep represents symmetry](../conventions.md) for the software's
    coordinate and array conventions. For executable examples, use the
    [tutorials](../examples/examples.md).

## Start with the representation you need

Begin with [space-group irreps](irreps/spacegroup_irreps.md) for reciprocal-space
conventions, little groups, and small representations. Then use the
[representation guide](irreps/irreps.md) to choose an extension for real vector
spaces, spinors, or magnetic symmetry. You can skip extensions that your calculation
does not use.

## Understand the calculation

The [algorithm guide](algorithm/algorithm.md) compares regular-representation and
subgroup-chain enumeration. Read it when choosing an algorithm or studying the
implementation. [Symmetry-adapted bases](projection.md) explains how projection
extracts basis functions for a chosen irrep.

## Further reading

For general background, see [Bradley and Cracknell](../references.md#Bradley2009-ze),
[Inui et al.](../references.md#Inui1996-et), and
[El-Batanouny and Wooten](../references.md#el2008symmetry).
Return to the [API reference](../api/api.md) for arguments and numerical settings.
