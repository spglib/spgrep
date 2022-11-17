# Change Log

## v0.3.0 (17 Nov. 2022)
- Add spinor representation [[#41]](https://github.com/spglib/spgrep/pull/41)

## v0.2.12
- Generate symmetric tensor by Erdös method {func}`spgrep.tensors.get_symmetry_adapted_tensors` and {func}`spgrep.tensors.apply_intrinsic_symmetry`
- Fix induced representation

## v0.2.11
- Add physically irreducible representation of space group
- Newly return Frobenius-Schur indicator from `enumerate_small_representations` and `enumerate_unitary_irreps`

## v0.2.9
- Add example notebooks

## v0.2.8

Initial release to PyPI
- Calculate irreducible representations (irreps) of space groups: {func}`spgrep.get_spacegroup_irreps` and {func}`spgrep.get_spacegroup_irreps_from_primitive_symmetry`
- Calculate irreps of crystallographic point groups: {func}`spgrep.get_crystallographic_pointgroup_irreps_from_symmetry`
- Calculate physically irreducible representations (irreps over real numbers)
- Apply projection operator: {func}`spgrep.representation.project_to_irrep`
- Unique decomposition of crystallographic point group: {func}`spgrep.pointgroup.get_pointgroup_chain_generators`
