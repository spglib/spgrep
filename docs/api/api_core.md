# Core functions

## Summary

Spgrep provides its core functions to enumerate irreducible representations (irreps) from given

- crystal structures (``(lattice, positions, numbers)``) or magnetic crystal structures (``(lattice, positions, numbers, magmoms)``),
- symmetry operations of space groups (``(rotations, translations)``) or magnetic space groups (``(rotations, translations, time_reversals)``),
- and symmetry operations of crystallographic point groups (``rotations``).

For a given crystal structure or symmetry operations, Spgrep calculates the following representations:

- Linear irreps (See [Irreps of space group](../formulation/irreps/spacegroup_irreps.md#space_group_irreps))
- Physically irreps (See [Physically irreducible representation (PIR)](../formulation/irreps/reality.md#physically_irreps))
- Projective irreps for spinor (See [Spin representation](../formulation/irreps/spinor.md#spin_representation))
- Projective irreducible co-representation (co-reps) for spinor (See [Co-representation of magnetic point and space group](../formulation/irreps/corep.md#corep))

The following tables summarize core functions for a combination of representations and inputs.

### Functions for (magnetic) crystal structure

| Tasks                         | (Magnetic) crystal structure                                |
|-------------------------------|------------------------------------------------------------ |
| Linear irreps                 | [`spgrep.get_spacegroup_irreps`][spgrep.get_spacegroup_irreps]                        |
| Projective irreps for spinor  | [`spgrep.get_spacegroup_spinor_irreps`][spgrep.get_spacegroup_spinor_irreps]                 |
| Projective co-reps for spinor | [`spgrep.get_spacegroup_spinor_irreps`][spgrep.get_spacegroup_spinor_irreps] with `magmoms`  |

### Functions for (magnetic) space-group symmetry operations

| Tasks                         | (Magnetic) space group                                                                    |
|-------------------------------|-------------------------------------------------------------------------------------------|
| Linear irreps                 | [`spgrep.get_spacegroup_irreps_from_primitive_symmetry`][spgrep.get_spacegroup_irreps_from_primitive_symmetry]                              |
| Physically irreps             | [`spgrep.get_spacegroup_irreps_from_primitive_symmetry`][spgrep.get_spacegroup_irreps_from_primitive_symmetry] with `real=True`             |
| Projective irreps for spinor  | [`spgrep.get_spacegroup_spinor_irreps_from_primitive_symmetry`][spgrep.get_spacegroup_spinor_irreps_from_primitive_symmetry]                       |
| Projective co-reps for spinor | [`spgrep.get_spacegroup_spinor_irreps_from_primitive_symmetry`][spgrep.get_spacegroup_spinor_irreps_from_primitive_symmetry] with `time_reversals` |

### Functions for (magnetic) crystallographic point-group symmetry operations

| Tasks                         | (Magnetic) point group                                                                           |
|-------------------------------|--------------------------------------------------------------------------------------------------|
| Linear irreps                 | [`spgrep.get_crystallographic_pointgroup_irreps_from_symmetry`][spgrep.get_crystallographic_pointgroup_irreps_from_symmetry]                              |
| Physically irreps             | [`spgrep.get_crystallographic_pointgroup_irreps_from_symmetry`][spgrep.get_crystallographic_pointgroup_irreps_from_symmetry] with `real=True`             |
| Projective irreps for spinor  | [`spgrep.get_crystallographic_pointgroup_spinor_irreps_from_symmetry`][spgrep.get_crystallographic_pointgroup_spinor_irreps_from_symmetry]                       |
| Projective co-reps for spinor | [`spgrep.get_crystallographic_pointgroup_spinor_irreps_from_symmetry`][spgrep.get_crystallographic_pointgroup_spinor_irreps_from_symmetry] with `time_reversals` |

## Linear representation

### Space group

::: spgrep.get_spacegroup_irreps

::: spgrep.get_spacegroup_irreps_from_primitive_symmetry

### Crystallographic point group

::: spgrep.get_crystallographic_pointgroup_irreps_from_symmetry

## Spin representation

### Space group

::: spgrep.get_spacegroup_spinor_irreps

::: spgrep.get_spacegroup_spinor_irreps_from_primitive_symmetry

### Crystallographic point group

::: spgrep.get_crystallographic_pointgroup_spinor_irreps_from_symmetry
