from spgrep.utils import is_matrix_group


def test_is_matrix_group(C3v):
    assert is_matrix_group(C3v)
