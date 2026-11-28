# Apply irreps to a physical problem

Choose a worked tutorial to construct symmetry-adapted bases or calculate spinor representations.

!!! info "Before you start"

    Complete the [quickstart](../quickstart.md) and read the
    [array conventions](../conventions.md). Each tutorial states its additional
    physics and Python prerequisites. For function signatures, use the
    [API reference](../api/api.md).

## Choose a tutorial

- [Construct symmetry-adapted tensors](symmetry_adapted_tensor.md) starts with cubic
  symmetry, decomposes a strain representation, and constructs invariant tensors.
- [Block-diagonalize a dynamical matrix](lattice_vibration.md) builds the displacement
  representation of a perovskite and projects a phonon dynamical matrix into irrep blocks.
- [Calculate spinor representations](spin_representation.md) uses the point group
  $C_{3v}$ to show how spin rotations and a factor system enter a projective representation.

Start with the tensor tutorial if you want to learn projection through a real-vector
example. The phonon and spinor tutorials can be read independently after the quickstart.

## Run a notebook

The pages display saved notebook outputs. Download links appear at the top of each
page. To run the tensor or spinor notebook, install the notebook and plotting tools:

```shell
python -m pip install spgrep notebook matplotlib seaborn
```

The phonon tutorial also needs Phonopy and the data file in the repository's
`examples` directory. Its prerequisites describe the checkout and working directory.
To understand the projection used in the tutorials, read
[symmetry-adapted bases](../formulation/projection.md).
