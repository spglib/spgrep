Formulation
===========

Working Example
---------------
Consider crystallographic point group :math:`3m` associated with :math:`P3m1` (No. 156).
The matrix representation of this crystallographic point group is

.. math::
    \mathcal{P}
    &=
    \left\{
        g_0=e,
        g_1,
        g_2=g_1^{-1},
        g_3,
        g_4=g_1^{-1} * g_3,
        g_5=g_1 * g_3
    \right\} \\
    \quad \mbox{where}\quad
    g_1
    &=
    \begin{pmatrix}
        0 & -1 & 0 \\
        1 & -1 & 0 \\
        0 & 0  & 1 \\
    \end{pmatrix} \\
    g_3
    &=
    \begin{pmatrix}
        0 & -1 & 0\\
        -1 & 0 & 0\\
        0 & 0 & 1\\
    \end{pmatrix}.

References
----------
.. [Net73] N. Neto, Acta Cryst. A, 29(4) 464–472 (1973).
.. [TVdV17] John C. Thomas and Anton Van der Ven, J. Mech. Phys. Solids 107, 76–95, (2017).
