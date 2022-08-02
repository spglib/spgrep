# Reality of (projective) irrep

Let {math}`(\Gamma, \mathrm{Span}_{\mathbb{C}} \{ \mathbf{v}_{i} \}_{i=1}^{d} )` be a unitary (projective) irrep of group {math}`G`.
Then, its conjugate (projective) representation is {math}`(\Gamma^{\ast}, \mathrm{Span}_{\mathbb{C}} \{ \mathbf{v}_{i}^{\ast} \}_{i=1}^{d} )`.
There are three cases for {math}`\Gamma` and {math}`\Gamma^{\ast}`:
1. {math}`\Gamma` is real: {math}`\Gamma` and {math}`\Gamma^{\ast}` are equivalent and can be taken as real matrices.
2. {math}`\Gamma` is pseudo-real: {math}`\Gamma` and {math}`\Gamma^{\ast}` are equivalent but can be taken as real matrices.
3. {math}`\Gamma` is not equivalent to {math}`\Gamma^{\ast}`.
These cases are classified with Frobenius-Schur indicator:
```{math}
  \frac{1}{|G|} \sum_{ g \in G } \chi(g^{2})
  = \begin{cases}
    1 & \mbox{($\Gamma$ is real)} \\
    -1 & \mbox{($\Gamma$ is pseudo-real)} \\
    0 & \mbox{($\Gamma$ is not equivalent to $\Gamma^{\ast}$)} \\
  \end{cases}.
```

We sometimes need to restrict irrep under a vector space over {math}`\mathbb{R}` (instead of {math}`\mathbb{C}`), which is called physically irreducible representation {cite}`PhysRevB.43.11010`.

## (1) {math}`\Gamma` is real

Refs. {cite}`Stokes:pc5025,Inui1996-et`

In this case, since {math}`\Gamma` and {math}`\Gamma^{\ast}` are equivalent, there exists a symmetric unitary matrix with
```{math}
  \mathbf{\Gamma}(g) \mathbf{U} &= \mathbf{U} \mathbf{\Gamma}(g)^{\ast} \\
  \mathbf{U}^{\dagger} \mathbf{U} &= \mathbf{1} \\
  \mathbf{U}^{\top} &= \mathbf{U}.
```
We impose {math}`\det \mathbf{U} = 1` additionally to specify the intertwiner uniquely.

The symmetric unitary matrix {math}`\mathbf{U}` can be diagonalized with real orthogonal matrix {math}`\mathbf{S}` as {math}`\mathbf{U} = \mathbf{S}^{-1} \mathbf{\Omega} \mathbf{S}` [^sym_unitary].
For this symmetric unitary matrix {math}`\mathbf{U}`, we can choose its square root with
```{math}
  \mathbf{T} &:= \mathbf{S}^{-1} \mathbf{\Omega}^{1/2} \mathbf{S} \\
  \mathbf{T}^{2} &= \mathbf{U} \\
  \mathbf{T}^{\dagger} \mathbf{T} &= \mathbf{1} \\
  \mathbf{T}^{\top} &= \mathbf{T}.
```

A transformed representation {math}`\mathbf{\Gamma}'(g) := \mathbf{T}\mathbf{\Gamma}(g)\mathbf{T}^{-1}` is real because
```{math}
  \mathbf{\Gamma}'(g)
  &= \mathbf{T} \mathbf{U}^{-1} \mathbf{\Gamma}(g)^{\ast} \mathbf{U} \mathbf{T}^{-1} \\
  &= \mathbf{T}^{-1} \mathbf{\Gamma}(g)^{\ast} \mathbf{T} \\
  &= \mathbf{\Gamma}'(g)^{\ast}.
```

[^sym_unitary]: For symmetric unitary matrix {math}`\mathbf{U}`, if {math}`\mathbf{v}` is eigenvector of {math}`\mathbf{U}`, its conjugacy {math}`\mathbf{v}'` is also eigenvector.
    Thus, we can take real part or imaginary part of {math}`\mathbf{v}` as a new basis vector.


## (2, 3) {math}`\Gamma` is pseudo-real or not equivalent to {math}`\Gamma^{\ast}`

transform conjugated basis pair to real vectors by unitary matrix:
```{math}
  (\mathbf{v}_{1}, \cdots, \mathbf{v}_{d}, \mathbf{v}_{1}^{\ast}, \cdots, \mathbf{v}_{d}^{\ast}) \mathbf{U}
    &= \sqrt{2} (\mathrm{Re}\, \mathbf{v}_{1}, \cdots, \mathrm{Re}\, \mathbf{v}_{d}, \mathrm{Im}\, \mathbf{v}_{1}, \cdots, \mathrm{Im}\, \mathbf{v}_{d}) \\
  \mathbf{U} &:= \frac{1}{\sqrt{2}}\begin{pmatrix}
    \mathbf{1}_{d} & -i \mathbf{1}_{d} \\
    \mathbf{1}_{d} & i \mathbf{1}_{d} \\
  \end{pmatrix} \quad (\mathrm{Unitary}) \\
  \mathbf{U}^{-1}
  \begin{pmatrix}
    \mathbf{D}(g) & \\
    & \mathbf{D}(g)^{\ast}
  \end{pmatrix}
  \mathbf{U}
  &= \begin{pmatrix}
    \mathrm{Re}\, \mathbf{D}(g) & \mathrm{Im}\, \mathbf{D}(g) \\
    -\mathrm{Im}\, \mathbf{D}(g) & \mathrm{Re}\, \mathbf{D}(g) \\
  \end{pmatrix}
  \quad (g \in G)
```

## References

```{bibliography}
:filter: docname in docnames
```
