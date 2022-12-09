"""Spin co-representation."""
from __future__ import annotations

import numpy as np

from spgrep.group import get_cayley_table
from spgrep.spinor import get_spinor_unitary_rotation
from spgrep.utils import NDArrayBool, NDArrayComplex, NDArrayFloat, NDArrayInt


def get_corep_spinor_factor_system(
    lattice: NDArrayFloat,
    rotations: NDArrayInt,
    time_reversals: NDArrayInt,
) -> tuple[NDArrayComplex, NDArrayComplex, NDArrayBool]:
    r"""Calculate spin-derived factor system of spin co-representation.

    See :ref:`corep_spinor_factor_system` for spinor-derived factor system :math:`\omega`.
    An ordinary symmetry operation :math:`\mathbf{S}_{i}1` maps to unitary operator :math:`\mathbf{U}(\mathbf{S}_{i})`.
    An antisymmetry operation :math:`\mathbf{S}_{i}1'` maps to

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    rotations: array[int], (order, 3, 3)
        Rotation parts of magnetic point group, :math:`\mathbf{S}_{i}`
    time_reversals: array[int], (order, )
        Time-reversal parts of magnetic point group, :math:`\theta_{i}`

    Returns
    -------
    corep_spinor_factor_system: array, (order, order)
        ``factor_system[i, j]`` stands for :math:`\omega(\mathbf{S}_{i}\theta_{i}, \mathbf{S}_{j}\theta_{j})`
    unitary_rotations: array, (order, 2, 2)
        ``unitary_rotations[i]`` stands for :math:`\mathbf{U}(\mathbf{S}_{i}) \in SU(2)`.
    anti_linear: array[bool], (order, )
        If ``anti_linear[i] == True``, the ``i``-th operator is anti-linear.
    """
    # -i sigma_y
    time_reversal_matrix = np.array(
        [
            [0, -1],
            [1, 0],
        ],
        dtype=np.complex_,
    )

    # Assign a unitary or anti-unitary operator for each magnetic operation
    order = len(rotations)
    unitary_rotations = np.zeros((order, 2, 2), dtype=np.complex_)
    anti_linear = np.zeros((order,), dtype=np.bool_)
    for i, (rotation, time_reversal) in enumerate(zip(rotations, time_reversals)):
        unitary_rotation = get_spinor_unitary_rotation(lattice, rotation)
        if time_reversal == 1:
            # Anti-linear operator
            unitary_rotations[i] = unitary_rotation @ time_reversal_matrix
            anti_linear[i] = True
        else:
            # Linear operator
            unitary_rotations[i] = unitary_rotation
            anti_linear[i] = False

    # Factor system for spinor co-rep
    table = get_cayley_table(rotations, time_reversals)
    corep_spinor_factor_system = np.zeros((order, order), dtype=np.complex_)
    for i, (ui, ai) in enumerate(zip(unitary_rotations, anti_linear)):
        for j, uj in enumerate(unitary_rotations):
            if ai:
                uiuj = ui @ np.conj(uj)  # ui is anti-linear
            else:
                uiuj = ui @ uj

            # si @ sj = sk in O(3)
            k = table[i, j]
            # Multiplier should be -1 or 1
            for multiplier in [-1, 1]:
                if np.allclose(uiuj, multiplier * unitary_rotations[k]):
                    corep_spinor_factor_system[i, j] = multiplier
                    break

    return corep_spinor_factor_system, unitary_rotations, anti_linear
