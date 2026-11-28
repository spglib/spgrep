# How spgrep represents symmetry

Understand how a crystal cell, a wave vector, and an ordered list of symmetry
operations determine the representation arrays returned by spgrep.

!!! info "Before you start"

    Work through the [quickstart](quickstart.md), or bring experience with NumPy
    arrays and crystallographic symmetry. This page explains the software's model.
    The [space-group theory](formulation/irreps/spacegroup_irreps.md) gives the derivation;
    the [API reference](api/api_core.md) specifies individual arguments and return values.

## A structure or a list of operations

A structure input consists of lattice vectors, fractional positions, and integer
species labels. spgrep uses spglib to find its symmetry operations. An operation
acts on a column of fractional coordinates as $x' = Rx + t$.
The arrays `rotations` and `translations` keep the pairs $(R,t)$ in the same order.

The structure-based functions accept primitive and conventional cells. Functions
with `from_primitive_symmetry` in their name instead require operations expressed
in a primitive cell. They do not search for symmetry or convert a conventional
cell on your behalf. Changing the cell basis also changes the coordinates of the
wave vector.

## A wave vector selects a little group

By default, `kpoint` gives coefficients in the reciprocal basis dual to the input
lattice. The lattice vectors are rows, and the reciprocal basis is
`np.linalg.inv(lattice).T`, without a factor of $2\pi$.
A lattice translation with fractional coordinates $n$ therefore has phase
$\exp(-2\pi i\,\mathrm{kpoint}\cdot n)$.

The little group contains the space-group operations that preserve the wave vector
modulo a reciprocal-lattice vector. spgrep returns matrices for finitely many
representatives of this group, with lattice translations accounted for by their phases.
These are the *small representations* used in the theory pages. The calculation
does not return an induced representation on every arm of the wave-vector star.

For centered input cells, the structure-based functions return operations and
matrices in the input cell, including its centering translations. Use the returned
mapping instead of assuming one operation per distinct rotation.

## Matrix indices follow operation indices

For the complex space-group calculation in the quickstart, each element of
`irreps` is an array with shape `(little_group_order, dim, dim)`.
The first axis indexes operations; the last two axes are matrix rows and columns.

`mapping_little_group[j]` gives the index in the full operation list corresponding
to matrix `irrep[j]`. Index `rotations` and `translations` with this mapping before
combining them with an irrep. The irrep list has no promised ordering, and a basis
change can alter matrix entries without changing the equivalence class.

A character is the trace of a representation matrix. Characters can help compare
equivalent representations after matching operations and their conventions.
An irrep dimension or its position in the list is not a crystallographic irrep label.

## Choose the representation space

Complex irreps
:   The default space-group and point-group calculations return unitary matrix
    representations on complex vector spaces.

Physically irreducible representations
:   Functions that accept `real=True` return representations on real vector spaces.
    This can combine complex irreps and increase the dimension. For primitive
    space-group input, the operation mapping then refers to the little group of
    $\pm\mathbf{k}$. See [physically irreducible representations](formulation/irreps/reality.md).

Spinor representations
:   Spinor functions return projective representations, a spinor factor system,
    and unitary spin rotations. The lattice is needed to convert spatial rotations
    to Cartesian coordinates. See [spin representation](formulation/irreps/spinor.md).

Co-representations
:   Magnetic symmetry can include anti-linear operations. A matrix and its
    anti-linearity flag must be interpreted together. See
    [co-representations](formulation/irreps/corep.md) and the relevant function's return values.

## Next steps

Use the [core API reference](api/api_core.md) to select a function for your input,
or [calculate irreps from primitive symmetry](howto/primitive-symmetry.md).
For symmetry-adapted bases, follow the [tensor tutorial](examples/symmetry_adapted_tensor.md).
