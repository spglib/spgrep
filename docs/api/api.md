# Find an API function

Locate the function reference for your input data and the calculation you need.

!!! info "Before you start"

    This reference assumes the [coordinate and array conventions](../conventions.md).
    If you have not run spgrep yet, start with the [quickstart](../quickstart.md).
    Function pages describe arguments and return values; the
    [theory guide](../formulation/formulation.md) contains the derivations.

## Start from a structure or symmetry operations

[Core functions](api_core.md) enumerate irreps from a crystal structure, primitive
space-group operations, or crystallographic point-group operations. Start there
for linear, real, or spinor representations.

## Work with representations you already have

- [Representation theory](api_rep.md) covers finite-group tables, characters,
  equivalence, decomposition, and enumeration from multiplication tables.
- [Crystallographic symmetry](api_symmetry.md) covers little groups, projection,
  symmetry-adapted tensors, cell transformations, and subgroup enumeration.
- [Spinor functions](api_spinor.md) construct spin rotations, spinor factor systems,
  and spinor small representations.
- [Co-representation functions](api_corep.md) handle spinor representations with
  anti-linear operations.
- [Utility functions](api_utils.md) includes the Grassmann distance for comparing subspaces.

## Next steps

Use the [primitive-symmetry how-to](../howto/primitive-symmetry.md) for a complete
operation-based calculation, or the [tutorials](../examples/examples.md) for applications.
