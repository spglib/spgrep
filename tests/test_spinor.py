import numpy as np
import pytest

from spgrep.group import check_cocycle_condition, get_cayley_table, get_identity_index
from spgrep.representation import is_representation, is_unitary
from spgrep.spinor import (
    enumerate_spinor_small_representations,
    get_rotation_angle_and_axis,
    get_spinor_factor_system_and_rotations,
)


@pytest.fixture
def hexagonal_lattice():
    a = 2.0
    c = 3.0
    lattice = np.array(
        [
            [a, 0, 0],
            [-0.5 * a, np.sqrt(3) / 2 * a, 0],
            [0, 0, c],
        ]
    )
    return lattice


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
    order = len(rotations)

    factor_system, unitary_rotations = get_spinor_factor_system_and_rotations(
        hexagonal_lattice,
        little_rotations=rotations,
        little_translations=np.zeros((order, 3)),
        kpoint=np.zeros(3),
    )

    # Check unitary
    for unitary_rotation in unitary_rotations:
        assert np.allclose(
            unitary_rotation @ np.conj(unitary_rotation).T,
            np.eye(2, dtype=np.complex_),
        )

    # Check factor system
    assert check_cocycle_condition(rotations, factor_system)
    table = get_cayley_table(rotations)
    identity_idx = get_identity_index(table)
    assert np.allclose(factor_system[identity_idx, :], 1)
    assert np.allclose(factor_system[:, identity_idx], 1)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
def test_spinor_irreps(method, C3v, hexagonal_lattice):
    # P3m1 (No. 156)
    rotations = C3v
    order = len(rotations)

    irreps, _ = enumerate_spinor_small_representations(
        hexagonal_lattice,
        little_rotations=rotations,
        little_translations=np.zeros((order, 3)),
        kpoint=np.zeros(3),
        method=method,
    )

    factor_system, _ = get_spinor_factor_system_and_rotations(
        hexagonal_lattice,
        little_rotations=rotations,
        little_translations=np.zeros((order, 3)),
        kpoint=np.zeros(3),
    )
    table = get_cayley_table(rotations)
    for irrep in irreps:
        assert is_representation(irrep, table, factor_system)
        assert is_unitary(irrep)

    # Check dimensions
    assert sorted([irrep.shape[1] for irrep in irreps]) == [1, 1, 2]
