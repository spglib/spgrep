import numpy as np

from spgrep.corep import get_corep_spinor_factor_system
from spgrep.group import get_cayley_table
from spgrep.utils import NDArrayComplex, NDArrayInt


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


def test_corep_spinor_factor_system(P42mnm_type3):
    rotations, _, time_reversals = P42mnm_type3
    lattice = np.eye(3)

    corep_spinor_factor_system, unitary_rotations, anti_linear = get_corep_spinor_factor_system(
        lattice, rotations, time_reversals
    )

    # Check unitary
    for unitary_rotation in unitary_rotations:
        assert np.allclose(
            unitary_rotation @ np.conj(unitary_rotation).T,
            np.eye(2, dtype=np.complex_),
        )

    # Cocycle condition
    assert check_corep_cocycle_condition(rotations, time_reversals, corep_spinor_factor_system)
