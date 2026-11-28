# Calculate your first space-group irreps

Calculate the irreducible representations of rutile at the X point and identify the
symmetry operations represented by each matrix.

!!! info "Before you start"

    You need Python and basic NumPy array notation, plus familiarity with a crystal
    cell and a reciprocal-space wave vector. Install spgrep using the
    [installation command](https://spglib.github.io/spgrep/#installation).
    Run the Python blocks below in order in one session.

This example uses complex representations without spin. For real representations,
spinors, and the mathematical derivation, follow the links at the end.

## 1. Describe the crystal

Import the calculation and character functions, then define the rutile cell:

```python
import numpy as np

from spgrep import get_spacegroup_irreps
from spgrep.rep.representation import get_character

# Rutile structure (https://materialsproject.org/materials/mp-2657/)
# P4_2/mnm (No. 136)
a = 4.603
c = 2.969
x_4f = 0.3046
lattice = [
    [a, 0, 0],
    [0, a, 0],
    [0, 0, c],
]
positions = [
    [0, 0, 0],  # Ti(2a)
    [0.5, 0.5, 0.5],  # Ti(2a)
    [x_4f, x_4f, 0],  # O(4f)
    [-x_4f, -x_4f, 0],  # O(4f)
    [-x_4f + 0.5, x_4f + 0.5, 0.5],  # O(4f)
    [x_4f + 0.5, -x_4f + 0.5, 0.5],  # O(4f)
]
numbers = [22, 22, 8, 8, 8, 8]
```

`lattice` contains lattice vectors as rows, `positions` contains fractional atomic
coordinates, and `numbers` labels atomic species. This cell is now ready for a
symmetry calculation.

## 2. Calculate irreps at X

Set the wave vector in the reciprocal basis of the input cell:

```python
kpoint = [0.5, 0, 0]
irreps, rotations, translations, mapping_little_group = get_spacegroup_irreps(
    lattice, positions, numbers, kpoint
)
print(len(rotations), len(mapping_little_group))
print(sorted(irrep.shape for irrep in irreps))
```

The output is:

```text
16 8
[(8, 2, 2), (8, 2, 2)]
```

The structure has 16 symmetry operations modulo translations of the input cell.
Eight belong to the little group at X, which preserves the wave vector modulo a
reciprocal-lattice vector. You now have two irreps, each with a 2-by-2 matrix for
each of those eight operations.

## 3. Match matrices to operations

Select the operations in the same order as the matrices:

```python
little_rotations = rotations[mapping_little_group]
little_translations = translations[mapping_little_group]

for irrep in irreps:
    assert irrep.shape == (len(little_rotations), 2, 2)
    character = get_character(irrep)
    identity = np.all(little_rotations == np.eye(3, dtype=int), axis=(1, 2))
    identity &= np.all(np.isclose(little_translations, 0), axis=1)
    print(character[identity].real)
```

Each line prints `[2.]`, the character of the identity for a two-dimensional irrep.
For any index `j`, `irrep[j]` represents
`(little_rotations[j], little_translations[j])`.
You can now associate the representation matrices and characters with crystal operations.

Equivalent irreps can have different bases and appear in a different order.
Compare characters for matching operations when checking results, rather than
expecting particular matrix entries or an irrep ordering.

## Next steps

Read [how spgrep represents symmetry](conventions.md) for coordinate and array conventions.
Use [primitive symmetry input](howto/primitive-symmetry.md) if you already have operations,
or choose a [worked tutorial](examples/examples.md) to apply the irreps.
For the derivation, start with [space-group irreps](formulation/irreps/spacegroup_irreps.md).
