from pathlib import Path

import numpy as np
import phonopy

from spgrep import get_spacegroup_irreps
from spgrep.representation import check_spacegroup_representation, project_to_irrep


def get_displacements_representation(
    lattice,
    positions,
    little_rotations,
    little_translations,
    qpoint,
):
    r"""Compute representation matrix for fourier-transformed displacements.

    .. math::
       \\Gamma_{\\kappa'\\mu'; \\kappa\\mu}^{\\mathbf{q}}(g) := \\exp \\left( -i \\mathbf{R}_{g} \\mathbf{q} \\cdot \\mathbf{h}_{g}(\\kappa) \\right) [\\mathbf{R}_{g}]_{\\mu'\\mu} \\delta_{ g\\kappa, \\kappa' }
    """
    little_order = len(little_rotations)
    num_atoms = len(positions)

    # Operation-`i` moves atom-`kappa` to `permutations[i, kappa]`
    permutations = np.zeros((little_order, num_atoms), dtype=int)
    for i, (Ri, vi) in enumerate(zip(little_rotations, little_translations)):
        for kappa, position in enumerate(positions):
            new_pos = np.remainder(Ri @ position + vi, 1)
            for kappa2, position2 in enumerate(positions):
                if np.allclose(position2, new_pos):
                    permutations[i, kappa] = kappa2
                    break

    shifts = np.zeros((little_order, num_atoms, 3))
    for i, (Ri, vi) in enumerate(zip(little_rotations, little_translations)):
        perm_i = permutations[i]
        shifts[i] = positions @ Ri.T + vi[None, :] - positions[perm_i]

    perm_rep = np.zeros((little_order, num_atoms, num_atoms), dtype=np.complex_)
    for i, Ri in enumerate(little_rotations):
        for kappa in range(num_atoms):
            kappa2 = permutations[i, kappa]
            perm_rep[i, kappa2, kappa] = np.exp(
                -2j * np.pi * np.dot(Ri.T @ qpoint, shifts[i, kappa])
            )

    # Rotation matrix in cartesian (order, 3, 3)
    A = np.transpose(lattice)  # column-wise lattice vectors
    Ainv = np.linalg.inv(A)
    rotation_rep = np.array([A @ r @ Ainv for r in little_rotations], dtype=np.complex_)

    rep = np.einsum("ipq,iab->ipaqb", perm_rep, rotation_rep, optimize="greedy")
    return rep.reshape(-1, num_atoms * 3, num_atoms * 3)


if __name__ == "__main__":
    # Perovskite structure: Pm-3m (No. 221)
    a = 3.986
    lattice = np.array(
        [
            [a, 0, 0],
            [0, a, 0],
            [0, 0, a],
        ]
    )
    positions = np.array(
        [
            [0, 0.5, 0.5],  # O(3c)
            [0.5, 0, 0.5],  # O(3c)
            [0.5, 0.5, 0],  # O(3c)
            [0.5, 0.5, 0.5],  # Ti(1b)
            [0, 0, 0],  # Ba(1a)
        ]
    )
    numbers = [0, 0, 0, 1, 2]

    qpoint = [0.5, 0, 0]  # X point (with primitive cell)
    irreps, rotations, translations, mapping_little_group = get_spacegroup_irreps(
        lattice, positions, numbers, qpoint
    )

    # Sanity check if `irreps` are representation for space group
    little_rotations = rotations[mapping_little_group]
    little_translations = translations[mapping_little_group]
    for irrep in irreps:
        assert check_spacegroup_representation(
            little_rotations, little_translations, qpoint, irrep
        )

    rep = get_displacements_representation(
        lattice, positions, little_rotations, little_translations, qpoint
    )
    assert check_spacegroup_representation(little_rotations, little_translations, qpoint, rep)

    all_basis = []
    for irrep in irreps:
        all_basis.extend(project_to_irrep(rep, irrep))

    path = Path(__file__).resolve().parent / "phonopy_mp-2998.yaml.xz"
    ph = phonopy.load(path)

    ph.dynamical_matrix.run(qpoint)
    dynamical_matrix = ph.dynamical_matrix.dynamical_matrix
