# Calculate irreps from primitive symmetry

Calculate complex small representations from space-group operations expressed in a primitive cell.

!!! info "Before you start"

    Install spgrep as in the [quickstart](../quickstart.md) and read the
    [coordinate conventions](../conventions.md).
    You need a complete set of operations modulo primitive lattice translations,
    with rotations, translations, and the wave vector in consistent bases.

For a crystal structure whose operations you have not yet found, use the
[structure-based quickstart](../quickstart.md). This procedure does not cover
spinors or magnetic symmetry.

## 1. Prepare operations in a primitive basis

Use your primitive-cell operations, or run this complete example for a group
containing identity and inversion:

```python
import numpy as np

from spgrep import get_spacegroup_irreps_from_primitive_symmetry

rotations = np.array([np.eye(3, dtype=int), -np.eye(3, dtype=int)])
translations = np.zeros((2, 3))
kpoint = np.array([0.5, 0, 0])
```

Each rotation acts on fractional coordinates. Each translation is expressed in
that same direct-lattice basis. `kpoint` uses the corresponding reciprocal basis.
If your operations describe a centered conventional cell, convert the cell,
operations, and wave vector consistently before calling this function.

## 2. Enumerate the representations

```python
irreps, mapping_little_group = get_spacegroup_irreps_from_primitive_symmetry(
    rotations, translations, kpoint
)
little_rotations = rotations[mapping_little_group]
little_translations = translations[mapping_little_group]
print(sorted(irrep.shape for irrep in irreps))
```

For this example the output is `[(2, 1, 1), (2, 1, 1)]`.

## 3. Use the returned operation mapping

To extract inversion characters, select inversion by its rotation and translation:

```python
inversion = np.all(little_rotations == -np.eye(3, dtype=int), axis=(1, 2))
inversion &= np.all(np.isclose(little_translations, 0), axis=1)
characters = [np.trace(irrep[inversion], axis1=1, axis2=2).item() for irrep in irreps]
print(sorted(round(value.real) for value in characters))
```

The output is `[-1, 1]`. When using another operation, select it in
`little_rotations` and `little_translations` rather than assuming its index.

## Next steps

Consult the [function reference](../api/api_core.md#spgrep.get_spacegroup_irreps_from_primitive_symmetry)
for tolerances and algorithm options. For real representations, read
[physically irreducible representations](../formulation/irreps/reality.md)
before setting `real=True`, which can change both the representation dimension and
the group selected by the mapping.
