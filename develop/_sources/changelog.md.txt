# Change Log

## v0.2.9
- Add example notebooks

## v0.2.8

Initial release to PyPI
- Calculate irreducible representations (irreps) of space groups: {func}`spgrep.get_spacegroup_irreps` and {func}`spgrep.get_spacegroup_irreps_from_primitive_symmetry`
- Calculate irreps of crystallographic point groups: {func}`spgrep.get_crystallographic_pointgroup_irreps_from_symmetry`
- Calculate physically irreducible representations (irreps over real numbers)
- Apply projection operator: {func}`spgrep.representation.project_to_irrep`
- Unique decomposition of crystallographic point group: {func}`spgrep.pointgroup.get_pointgroup_chain_generators`
