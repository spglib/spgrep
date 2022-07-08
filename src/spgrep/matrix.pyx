# Ref: https://cython.readthedocs.io/en/latest/src/tutorial/numpy.html
import numpy as np

cimport cython
cimport numpy as np

np.import_array()

@cython.boundscheck(False)
@cython.wraparound(False)
def block_diagonalize(
    np.complex128_t [:, :] reg_nonzero,
    int [:, :] lookup,
    np.complex128_t [:, :] transformation,
):
    # reg_nonzero: array, (order, order)
    #   reg_nonzero[m, i] == reg[m, i, lookup[m, i]]
    # lookup: array, (order, order)
    #   For (m, i), reg[m, i, lookup[m, i]] has only one nonzero entry.
    # transformation: array, (order, dim)
    cdef size_t i, j, k, l, m
    cdef double complex t_li_conj_reg_klm, reg_klm

    cdef size_t order = transformation.shape[0]
    cdef size_t dim = transformation.shape[1]
    cdef np.ndarray irrep = np.zeros([order, dim, dim], dtype=np.complex128)

    # irrep = np.einsum("li,klm,mj->kij", np.conj(transformation), reg, transformation)
    for k in range(order):
        for l in range(order):
            m = lookup[k, l]
            reg_klm = reg_nonzero[k, l]
            for i in range(dim):
                t_li_conj_reg_klm = np.conj(transformation[l, i]) * reg_klm
                for j in range(dim):
                    irrep[k, i, j] += t_li_conj_reg_klm * transformation[m, j]

    return irrep
