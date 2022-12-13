import numpy as np
import pytest

from spgrep.core import (
    get_crystallographic_pointgroup_spinor_irreps_from_symmetry,
    get_spacegroup_spinor_irreps,
    get_spacegroup_spinor_irreps_from_primitive_symmetry,
)
from spgrep.group import (
    check_cocycle_condition,
    get_cayley_table,
    get_identity_index,
    get_little_group,
)
from spgrep.representation import (
    check_spacegroup_representation,
    is_representation,
    is_unitary,
)
from spgrep.spinor import (
    enumerate_spinor_small_representations,
    get_rotation_angle_and_axis,
    get_spinor_factor_system,
)


@pytest.mark.parametrize(
    "cart_rotation,angle,cart_axis",
    [
        (np.diag([-1, -1, 1]), np.pi, np.array([0, 0, 1])),
        (np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]), np.pi / 2, np.array([0, 0, 1])),
    ],
)
def test_rotation_angle_and_axis(cart_rotation, angle, cart_axis):
    angle_actual, cart_axis_actual = get_rotation_angle_and_axis(cart_rotation)
    assert np.isclose(angle_actual, angle)
    assert np.allclose(cart_axis_actual, cart_axis)


def test_spinor_factor_system_symmorphic(C3v, hexagonal_lattice):
    # P3m1 (No. 156)
    rotations = C3v
    lattice = hexagonal_lattice

    spinor_factor_system, unitary_rotations = get_spinor_factor_system(lattice, rotations)

    # Check unitary
    for unitary_rotation in unitary_rotations:
        assert np.allclose(
            unitary_rotation @ np.conj(unitary_rotation).T,
            np.eye(2, dtype=np.complex_),
        )

    # Check factor system
    assert check_cocycle_condition(rotations, spinor_factor_system)
    table = get_cayley_table(rotations)
    identity_idx = get_identity_index(table)
    assert np.allclose(spinor_factor_system[identity_idx, :], 1)
    assert np.allclose(spinor_factor_system[:, identity_idx], 1)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
def test_spinor_irreps(method, C3v, hexagonal_lattice):
    # P3m1 (No. 156)
    rotations = C3v
    lattice = hexagonal_lattice

    irreps, spinor_factor_system, _ = enumerate_spinor_small_representations(
        lattice=lattice,
        little_rotations=rotations,
        method=method,
    )

    table = get_cayley_table(rotations)
    for irrep in irreps:
        # Factor system from nonsymmorphic is trivial in this case
        assert is_representation(irrep, table, spinor_factor_system)
        assert is_unitary(irrep)

    # Check dimensions
    assert sorted([irrep.shape[1] for irrep in irreps]) == [1, 1, 2]


@pytest.mark.parametrize(
    "kpoint,shape_expect",
    [
        ([0, 1, 1 / 2], [2, 2, 2]),  # T point for hR
        ([-1 / 2, 1 / 2, 1 / 2], [2]),  # L point for hR
    ],
)
def test_get_spacegroup_spinor_irreps(kpoint, shape_expect, corundum_cell):
    (
        irreps,
        little_spinor_factor_system,
        little_unitary_rotations,
        rotations,
        translations,
        mapping,
    ) = get_spacegroup_spinor_irreps(
        *corundum_cell,
        kpoint=kpoint,
    )
    assert [irrep.shape[1] for irrep in irreps] == shape_expect


@pytest.mark.parametrize(
    "kpoint,shape_expect",
    [
        ([0, 0, 0], [2, 2, 2, 2]),  # Gammma point
        ([0, 1 / 2, 0], [2, 2]),  # X point
        ([0, 0, 1 / 2], [4]),  # Z point
    ],
)
def test_get_spacegroup_spinor_irreps_from_primitive_symmetry(kpoint, shape_expect, P42mnm):
    rotations, translations = P42mnm
    (
        irreps,
        little_spinor_factor_system,
        little_unitary_rotations,
        mapping_little_group,
    ) = get_spacegroup_spinor_irreps_from_primitive_symmetry(
        lattice=np.eye(3),
        rotations=rotations,
        translations=translations,
        kpoint=kpoint,
    )

    # Check unitary rotations
    for unitary_rotation in little_unitary_rotations:
        assert np.allclose(
            unitary_rotation @ np.conj(unitary_rotation).T,
            np.eye(2, dtype=np.complex_),
        )

    # Check as representation
    little_rotations, little_translations, _ = get_little_group(rotations, translations, kpoint)
    for irrep in irreps:
        assert check_spacegroup_representation(
            little_rotations, little_translations, kpoint, irrep, little_spinor_factor_system
        )

    # Check dimensions of irreps
    assert [irrep.shape[1] for irrep in irreps] == shape_expect


def test_get_crystallographic_pointgroup_spinor_irreps_from_symmetry(Oh):
    rotations = Oh

    (
        irreps,
        factor_system,
        unitary_rotations,
    ) = get_crystallographic_pointgroup_spinor_irreps_from_symmetry(
        lattice=np.eye(3),
        rotations=rotations,
    )

    # Check unitary
    for unitary_rotation in unitary_rotations:
        assert np.allclose(
            unitary_rotation @ np.conj(unitary_rotation).T,
            np.eye(2, dtype=np.complex_),
        )

    # Check as representation
    table = get_cayley_table(rotations)
    for irrep in irreps:
        assert is_representation(irrep, table, factor_system)
        assert is_unitary(irrep)

    # Check irreps
    # Ref: https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_out.pl?tipogrupo=dbg&pointspace=point&num=221&super=32&symbol=m-3m
    assert sorted([irrep.shape[1] for irrep in irreps]) == [2, 2, 2, 2, 4, 4]
