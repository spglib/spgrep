import numpy as np
import pytest

from spgrep.core import (
    get_crystallographic_pointgroup_spinor_irreps_from_symmetry,
    get_magnetic_spacegroup_coreps_from_primitive_symmetry,
)
from spgrep.corep import get_corep_spinor_factor_system
from spgrep.group import get_cayley_table
from spgrep.utils import NDArrayComplex, NDArrayFloat, NDArrayInt


def check_corep_cocycle_condition(
    rotations: NDArrayInt,
    time_reversals: NDArrayInt,
    corep_factor_system: NDArrayComplex,
) -> bool:
    """Return true if given factor system satisfies the cocycle condition for co-representation."""
    order = len(rotations)
    table = get_cayley_table(rotations, time_reversals)

    # Check associativity: (gi * gj) * gk = gi * (gj * gk)
    for i, tri in enumerate(time_reversals):
        for j in range(order):
            for k in range(order):
                ij = table[i, j]
                jk = table[j, k]

                if tri == 1:
                    if not np.isclose(
                        corep_factor_system[i, j] * corep_factor_system[ij, k],
                        corep_factor_system[i, jk] * np.conj(corep_factor_system[j, k]),
                    ):
                        return False
                else:
                    if not np.isclose(
                        corep_factor_system[i, j] * corep_factor_system[ij, k],
                        corep_factor_system[i, jk] * corep_factor_system[j, k],
                    ):
                        return False

    return True


@pytest.mark.parametrize(
    "symmetry_and_lattice",
    [
        ("P42mnm_type3"),
        # ("bcc_type4"),  # Passed but too long to test
    ],
)
def test_corep_spinor_factor_system(request, symmetry_and_lattice):
    rotations, _, time_reversals, lattice = request.getfixturevalue(symmetry_and_lattice)

    corep_spinor_factor_system, unitary_rotations, anti_linear = get_corep_spinor_factor_system(
        lattice, rotations, time_reversals
    )

    # Check unitary
    for unitary_rotation in unitary_rotations:
        assert np.allclose(
            unitary_rotation @ np.conj(unitary_rotation).T,
            np.eye(2, dtype=np.complex128),
        )

    # Cocycle condition
    assert check_corep_cocycle_condition(rotations, time_reversals, corep_spinor_factor_system)


@pytest.mark.parametrize("method", [("Neto"), ("random")])
@pytest.mark.parametrize(
    "symmetry_and_lattice",
    [
        ("P42mnm_type1"),
        ("P42mnm_type2"),
        ("P42mnm_type3"),
        ("bcc_type4"),
    ],
)
def test_get_crystallographic_pointgroup_spinor_irreps_from_symmetry(
    request, method, symmetry_and_lattice
):
    rotations, _, time_reversals, lattice = request.getfixturevalue(symmetry_and_lattice)

    # TODO: Add more tests
    (
        co_irreps,
        indicators,
        factor_system,
        unitary_rotations,
        anti_linear,
    ) = get_crystallographic_pointgroup_spinor_irreps_from_symmetry(
        lattice,
        rotations,
        time_reversals,
        method=method,
    )


@pytest.mark.parametrize(
    "kpoint",
    [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.0, 0.0]),
        np.array([0.5, 0.5, 0.0]),
        np.array([0.5, 0.5, 0.5]),
    ],
)
@pytest.mark.parametrize("method", [("Neto"), ("random")])
@pytest.mark.parametrize(
    "symmetry_and_lattice",
    [
        ("P42mnm_type1"),
        ("P42mnm_type2"),
        ("P42mnm_type3"),
        ("bcc_type4"),
    ],
)
def test_get_magnetic_spacegroup_coreps_from_primitive_symmetry(
    request, kpoint: NDArrayFloat, method, symmetry_and_lattice
):
    rotations, translations, time_reversals, lattice = request.getfixturevalue(
        symmetry_and_lattice
    )

    # TODO: Add more tests
    (
        co_irreps,
        anti_linear,
        mapping_little_group,
    ) = get_magnetic_spacegroup_coreps_from_primitive_symmetry(
        rotations=rotations,
        translations=translations,
        time_reversals=time_reversals,
        kpoint=kpoint,
        method=method,
    )
