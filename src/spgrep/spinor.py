"""Spinor representation."""
from __future__ import annotations

from itertools import product
from typing import Literal

import numpy as np

from spgrep.group import get_cayley_table, get_factor_system_from_little_group
from spgrep.irreps import enumerate_unitary_irreps
from spgrep.utils import NDArrayComplex, NDArrayFloat, NDArrayInt


def enumerate_spinor_small_representations(
    lattice: NDArrayFloat,
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
    method: Literal["Neto", "random"] = "Neto",
    rtol: float = 1e-5,
    atol: float = 1e-8,
    max_num_random_generations: int = 4,
) -> tuple[list[NDArrayComplex], list[NDArrayComplex]]:
    """Enumerate all unitary irreps of little group for spinor.

    See :ref:`spinor_factor_system` for Spgrep's convention of spin-derived factor system :math:`z(S_{i}, S_{j})`.

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    little_rotations: array, (order, 3, 3)
    little_translations: array, (order, 3)
    kpoint: array, (3, )
    method: str, 'Neto' or 'random'
        'Neto': construct irreps from a fixed chain of subgroups of little co-group
        'random': construct irreps by numerically diagonalizing a random matrix commute with regular representation
    rtol: float
        Relative tolerance to distinguish difference eigenvalues
    atol: float
        Relative tolerance to compare
    max_num_random_generations: int
        Maximum number of trials to generate random matrix

    Returns
    -------
    irreps: list of unitary small representations (irreps of little group) with (order, dim, dim)
    unitary_rotations: array, (order, 2, 2)
        SU(2) rotations on spinor.
    """
    factor_system, unitary_rotations = get_spinor_factor_system_and_rotations(
        lattice, little_rotations, little_translations, kpoint
    )

    # Compute irreps of little co-group
    little_cogroup_irreps, _ = enumerate_unitary_irreps(
        little_rotations,
        factor_system,
        real=False,  # Nonsense to consider real-value irreps
        method=method,
        rtol=rtol,
        atol=atol,
        max_num_random_generations=max_num_random_generations,
    )

    # Small representations of little group
    irreps = []
    for rep in little_cogroup_irreps:
        phases = np.array(
            [
                np.exp(-2j * np.pi * np.dot(kpoint, translation))
                for translation in little_translations
            ]
        )
        irreps.append(rep * phases[:, None, None])

    return irreps, unitary_rotations


def get_spinor_factor_system_and_rotations(
    lattice: NDArrayFloat,
    little_rotations: NDArrayInt,
    little_translations: NDArrayFloat,
    kpoint: NDArrayFloat,
) -> tuple[NDArrayComplex, NDArrayComplex]:
    r"""Calculate factor system of spinor representation for little co-group.

    .. math::
       D^{\mathbf{k}\alpha}(S_{i}) D^{\mathbf{k}\alpha}(S_{j})
       = z(S_{i}, S_{j}) \exp \left( -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} \right) D^{\mathbf{k}\alpha}(S_{k})

    where :math:`S_{i}S_{j} = S_{k}` and :math:`\mathbf{g}_{i} = S_{i}^{-1} \mathbf{k} - \mathbf{k}`.
    See :ref:`spinor_factor_system` for Spgrep's convention of spin-derived factor system :math:`z(S_{i}, S_{j})`.

    Parameters
    ----------
    lattice: array, (3, 3)
        Row-wise basis vectors. ``lattice[i, :]`` is the i-th lattice vector.
    little_rotations: array, (order, 3, 3)
        Linear parts of coset of little group stabilizing ``kpoint``.
    little_translations: array, (order, 3)
        Translation parts of coset of little group stabilizing ``kpoint``.
    kpoint: array, (3, )

    Returns
    -------
    factor_system: array, (order, order)
        Factor system of spinor representations for little co-group that have one-to-one correspondence to small representations
    unitary_rotations: array, (order, 2, 2)
        SU(2) rotations on spinor.
    """
    order = len(little_rotations)

    # Assign a SU(2) rotation for each O(3) operation
    unitary_rotations = np.zeros((order, 2, 2), dtype=np.complex_)
    for i, si in enumerate(little_rotations):
        # Ignore inversion parts
        if np.linalg.det(si) > 0:
            sign = 1
        else:
            sign = -1
        rotation = si * sign

        cart_rotation = lattice.T @ rotation @ np.linalg.inv(lattice.T)
        theta, cart_axis = get_rotation_angle_and_axis(cart_rotation)
        cos_half = np.cos(0.5 * theta)
        sin_half = np.sin(0.5 * theta)
        unitary_rotations[i] = np.array(
            [
                [
                    cos_half - 1j * cart_axis[2] * sin_half,
                    sin_half * (-1j * cart_axis[0] - cart_axis[1]),
                ],
                [
                    sin_half * (-1j * cart_axis[0] + cart_axis[1]),
                    cos_half + 1j * cart_axis[2] * sin_half,
                ],
            ]
        )

    # Factor system from spin
    table = get_cayley_table(little_rotations)
    spin_factor_system = np.zeros((order, order))
    for (i, si), (j, sj) in product(enumerate(little_rotations), repeat=2):
        # si @ sj = sk in O(3)
        k = table[i, j]
        uiuj = unitary_rotations[i] @ unitary_rotations[j]
        # Multiplier should be -1 or 1
        for multiplier in [-1, 1]:
            if np.allclose(uiuj * multiplier, unitary_rotations[k]):
                spin_factor_system[i, j] = multiplier
                break

    # Factor system from nonsymmorphic
    nonsymmorphic_factor_system = get_factor_system_from_little_group(
        little_rotations=little_rotations,
        little_translations=little_translations,
        kpoint=kpoint,
    )

    factor_system = spin_factor_system * nonsymmorphic_factor_system
    return factor_system, unitary_rotations


def get_rotation_angle_and_axis(cart_rotation: NDArrayFloat) -> tuple[float, NDArrayFloat]:
    """Return angle and axis of a rotation.

    Angle is chosen between 0 and pi.
    """
    identity = np.eye(3)
    if np.allclose(cart_rotation, identity):
        # We choose theta=0 for identity operation as convention
        return 0.0, np.array([0, 0, 1])

    # Here, 0 < theta <= pi
    cos_theta = (np.trace(cart_rotation) - 1) / 2
    if np.isclose(cos_theta, -1):
        theta = np.pi
        eigvals, eigvecs = np.linalg.eig(cart_rotation)
        cart_axis = None
        for i in range(3):
            if not np.isclose(eigvals[i], 1):
                continue
            cart_axis = np.real(eigvecs[:, i])
            cart_axis /= np.linalg.norm(cart_axis)

            # Fix direction by lexicographic order
            for j in range(3):
                if not np.isclose(cart_axis[j], 0):
                    if cart_axis[j] < 0:
                        cart_axis *= -1
                    break
    else:
        nondiag = np.array(
            [
                cart_rotation[1, 2] - cart_rotation[2, 1],
                cart_rotation[2, 0] - cart_rotation[0, 2],
                cart_rotation[0, 1] - cart_rotation[1, 0],
            ]
        )
        sin_theta = np.linalg.norm(nondiag) / 2
        # Since sin_theta >= 0, theta in (0, pi)
        theta = np.arctan2(sin_theta, cos_theta)

        cart_axis = nondiag / (-2 * sin_theta)

    return theta, cart_axis
