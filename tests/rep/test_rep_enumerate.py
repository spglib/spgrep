from __future__ import annotations

from itertools import permutations

import numpy as np
import pytest

from spgrep.rep.enumerate import (
    _get_irreps_from_matrix,
    decompose_representation,
    enumerate_unitary_irreps_from_regular_representation,
)
from spgrep.rep.representation import get_character, is_representation, is_unitary


@pytest.mark.parametrize(
    "eigenvalues, rtol",
    [
        ([1.0, 1.075, 1.15], 0.1),
        ([-1.0, -0.925, -0.85], 0.1),
        ([0.0, 0.75e-8, 1.5e-8], 0.0),
    ],
)
def test_eigenspace_tolerance_uses_first_eigenvalue(eigenvalues, rtol):
    # Adjacent eigenvalues are close, but the first and last are not. The first
    # two trivial copies form a reducible space; the sign irrep must survive.
    representation = np.array([np.eye(3), np.diag([1, 1, -1])], dtype=complex)
    irreps = _get_irreps_from_matrix(representation, np.diag(eigenvalues), rtol=rtol)

    assert len(irreps) == 1
    for irrep in irreps:
        np.testing.assert_allclose(get_character(irrep), [1, -1])


@pytest.mark.parametrize("layout", ["C", "F", "strided"])
def test_decompose_dense_non_regular_representation(rng, layout):
    # The permutation representation of S3 plus its sign irrep contains each
    # of the three irreps once, in a four-dimensional space for a group of order six.
    representation = np.zeros((6, 4, 4), dtype=complex)
    for i, permutation in enumerate(permutations(range(3))):
        representation[i, :3, :3] = np.eye(3)[:, permutation]
        representation[i, 3, 3] = np.linalg.det(representation[i, :3, :3])
    basis, _ = np.linalg.qr(rng.random((4, 4)) + 1j * rng.random((4, 4)))
    representation = basis.conj().T @ representation @ basis
    if layout == "strided":
        padded = np.zeros((6, 8, 8), dtype=complex)
        padded[:, ::2, ::2] = representation
        representation = padded[:, ::2, ::2]
    else:
        representation = np.array(representation, order=layout)

    irreps = decompose_representation(representation)

    assert [irrep.shape[1] for irrep in irreps] == [1, 1, 2]
    np.testing.assert_allclose(
        [get_character(irrep) for irrep in irreps],
        [[1, 1, 1, 1, 1, 1], [1, -1, -1, 1, 1, -1], [2, 0, 0, -1, -1, 0]],
        atol=1e-8,
    )
    assert all(is_unitary(irrep) for irrep in irreps)


@pytest.mark.parametrize("layout", ["C", "F", "strided"])
def test_complex_projective_regular_representation(layout):
    # C2 x C2 with a Pauli cocycle has one two-dimensional projective irrep.
    elements = np.arange(4)
    table = elements[:, None] ^ elements[None, :]
    phases = np.exp(1j * np.array([0.0, 0.3, 0.7, 1.1]))
    factor_system = (
        (-1.0) ** ((elements[:, None] // 2) * (elements[None, :] % 2))
        * phases[:, None]
        * phases[None, :]
        / phases[table]
    )
    representation = np.zeros((4, 4, 4), dtype=complex)
    for i in elements:
        representation[i, table[i], elements] = factor_system[i]
    if layout == "strided":
        padded = np.zeros((4, 8, 8), dtype=complex)
        padded[:, ::2, ::2] = representation
        representation = padded[:, ::2, ::2]
    else:
        representation = np.array(representation, order=layout)

    irreps = enumerate_unitary_irreps_from_regular_representation(representation)

    assert [irrep.shape[1] for irrep in irreps] == [2]
    for irrep in irreps:
        np.testing.assert_allclose(get_character(irrep), [2, 0, 0, 0], atol=1e-8)
        assert is_unitary(irrep)
        assert is_representation(irrep, table, factor_system)


@pytest.mark.parametrize(
    "matrix",
    [
        pytest.param([[1, 0], [0, 0]], id="zero-row"),
        pytest.param([[1, 1], [0, 1]], id="multiple-entries"),
        pytest.param([[1, 1], [0, 0]], id="redistributed-nonzeros"),
    ],
)
def test_regular_representation_requires_one_nonzero_per_row(matrix):
    representation = np.array([matrix, [[0, 1], [1, 0]]], dtype=complex)

    with pytest.raises(ValueError):
        enumerate_unitary_irreps_from_regular_representation(representation)


def test_complex_conjugate_irreps_are_distinct():
    elements = np.arange(4)
    representation = np.zeros((4, 4, 4), dtype=int)
    for i in elements:
        representation[i, (i + elements) % 4, elements] = 1
    irreps = enumerate_unitary_irreps_from_regular_representation(representation)

    assert [irrep.shape[1] for irrep in irreps] == [1, 1, 1, 1]
    characters = np.array([get_character(irrep) for irrep in irreps])
    for expected in [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1j, -1, -1j], [1, -1j, -1, 1j]]:
        matches = np.all(np.isclose(characters, expected), axis=1)
        assert np.count_nonzero(matches) == 1
