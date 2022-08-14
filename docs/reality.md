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

We sometimes need to restrict irrep under a vector space over {math}`\mathbb{R}` (instead of {math}`\mathbb{C}`), which is called physically irreducible representation (PIR) {cite}`PhysRevB.43.11010`.

## PIR of finite group

### (1) {math}`\Gamma` is real

Refs. {cite}`Inui1996-et`

In this case, since {math}`\Gamma` and {math}`\Gamma^{\ast}` are equivalent, there exists a symmetric unitary matrix with
```{math}
  \mathbf{\Gamma}(g) \mathbf{U} &= \mathbf{U} \mathbf{\Gamma}(g)^{\ast} \\
  \mathbf{U}^{\dagger} \mathbf{U} &= \mathbf{1} \\
  \mathbf{U}^{\top} &= \mathbf{U}.
```
We impose {math}`\det \mathbf{U} = 1` additionally to specify the intertwiner uniquely.
The intertwiner can be numerically computed as shown in {ref}`here <intertwiner>`.

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

### (2, 3) {math}`\Gamma` is pseudo-real or not equivalent to {math}`\Gamma^{\ast}`

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

## PIR of space group {math}`\mathcal{G}`

Ref. {cite}`Stokes:pc5025`

Consider to construct PIR of space group from small representation {math}`\Gamma^{\mathbf{k}\alpha}` at {math}`\mathbf{k}`.

### Frobenius-Schur indicator for space-group representations

Let {math}`\Gamma^{(\mathbf{k}, 0)}` and {math}`\Gamma^{(\mathbf{k}, 1)}` be representations of space group {math}`\mathcal{G}` with {math}`\mathbf{k}` vector.
Let {math}`\mathcal{T}` be a translational subgroup of {math}`\mathcal{G}`.
Consider coset representatives of little group of {math}`\mathbf{k}` over {math}`\mathcal{T}`:
```{math}
  \mathcal{G}^{\mathbf{k}} = \coprod_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } (\mathbf{S}_{i}, \mathbf{w}_{i}) \mathcal{T}.
```

Then sum over translations in {math}`\mathcal{G}^{\mathbf{k}}`:
```{math}
  \frac{1}{|G|} \sum_{ g \in G } \chi^{\mathbf{k}}(g^{2})
  &= \frac{1}{N} \sum_{ \mathbf{t} } \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
      \chi^{\mathbf{k}}\left( (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) \right) \\
  &= \frac{1}{N} \sum_{ \mathbf{t} } \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
      \chi^{\mathbf{k}}\left( (\mathbf{E}, \mathbf{t} + \mathbf{S}_{i}\mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i})^{2} \right) \\
  &= \left(
        \frac{1}{N} \sum_{ \mathbf{t} } e^{-i \mathbf{k} \cdot (\mathbf{t} + \mathbf{S}_{i}\mathbf{t}) }
     \right)
     \left(
        \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } \chi^{\mathbf{k}}\left( (\mathbf{S}_{i}, \mathbf{w}_{i})^{2} \right) \\
     \right) \\
  &= \mathbb{I}\left[ 2\mathbf{k} \equiv \mathbf{0} \right] \cdot
      \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } \chi^{\mathbf{k}}\left( (\mathbf{S}_{i}, \mathbf{w}_{i})^{2} \right),
```
where {math}`\mathbb{I}[C]` takes one if the condition {math}`C` is true and takes zero otherwise.
Here we use {math}`\mathbf{S}_{i}^{\top} \in \overline{\mathcal{G}}^{\mathbf{k}}` as
```{math}
  \frac{1}{N} \sum_{ \mathbf{t} } e^{-i \mathbf{k} \cdot (\mathbf{t} + \mathbf{S}_{i}\mathbf{t}) }
    &= \frac{1}{N} \sum_{ \mathbf{t} } e^{-i (\mathbf{E} + \mathbf{S}_{i})^{\top} \mathbf{k} \cdot \mathbf{t} } \\
    &= \mathbb{I} \left[ (\mathbf{E} + \mathbf{S}_{i})^{\top} \mathbf{k} \equiv \mathbf{0} \right] \\
    &= \mathbb{I} \left[ \mathbf{k} + \mathbf{S}_{i}^{\top} \mathbf{k} \equiv \mathbf{0} \right] \\
    &= \mathbb{I} \left[ 2\mathbf{k} \equiv \mathbf{0} \right].
```

### (1) {math}`\Gamma^{\mathbf{k}\alpha}` is real

We need to find an intertwiner
```{math}
  \mathbf{\Gamma}^{\mathbf{k}\alpha}(g) \mathbf{U} = \mathbf{U} \mathbf{\Gamma}^{\mathbf{k}\alpha}(g)^{\ast}
```

Because we can assume {math}`2\mathbf{k} \equiv \mathbf{0}` in this case, the following matrix is an intertwiner:
```{math}
    \mathbf{U}
    &:= \frac{1}{N} \sum_{ \mathbf{t} } \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) \right)
            \mathbf{B}
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) \right)^{\ast \dagger} \\
    &= \left(
            \frac{1}{N} \sum_{ \mathbf{t} } e^{2i\mathbf{k}\cdot\mathbf{t}}
       \right)
       \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)
            \mathbf{B}
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)^{\ast \dagger} \\
    &= \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)
            \mathbf{B}
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)^{\ast \dagger}
       \quad (\because 2\mathbf{k} \cdot \mathbf{t} \in \mathbb{Z}) \\
```


### (2, 3) {math}`\Gamma^{\mathbf{k}\alpha}` is pseudo-real or not equivalent to {math}`\Gamma^{\mathbf{k}\alpha \ast}`

transform conjugated basis pair to real vectors by unitary matrix:
```{math}
  (\mathbf{v}_{1}, \cdots, \mathbf{v}_{d}, \mathbf{v}_{1}^{\ast}, \cdots, \mathbf{v}_{d}^{\ast}) \mathbf{U}
    &= \sqrt{2} (\mathrm{Re}\, \mathbf{v}_{1}, \cdots, \mathrm{Re}\, \mathbf{v}_{d}, \mathrm{Im}\, \mathbf{v}_{1}, \cdots, \mathrm{Im}\, \mathbf{v}_{d}) \\
  \mathbf{U} &:= \frac{1}{\sqrt{2}}\begin{pmatrix}
    \mathbf{1}_{d} & -i \mathbf{1}_{d} \\
    \mathbf{1}_{d} & i \mathbf{1}_{d} \\
  \end{pmatrix} \quad (\mathrm{Unitary}) \\
  \tilde{\mathbf{\Gamma}}^{\mathbf{k}\alpha}(g)
  &:=
  \mathbf{U}^{-1}
  \begin{pmatrix}
    \mathbf{\Gamma}^{\mathbf{k}\alpha}(g) & \\
    & \mathbf{\Gamma}^{\mathbf{k}\alpha}(g)^{\ast}
  \end{pmatrix}
  \mathbf{U} \\
  &= \begin{pmatrix}
    \mathrm{Re}\, \mathbf{\Gamma}^{\mathbf{k}\alpha}(g) & \mathrm{Im}\, \mathbf{\Gamma}^{\mathbf{k}\alpha}(g) \\
    -\mathrm{Im}\, \mathbf{\Gamma}^{\mathbf{k}\alpha}(g) & \mathrm{Re}\, \mathbf{\Gamma}^{\mathbf{k}\alpha}(g) \\
  \end{pmatrix}
  \quad (g \in \mathcal{G})
```

In particular, for translation {math}`(E, \mathbf{t}) \in \mathcal{G}`,
```{math}
  \tilde{\mathbf{\Gamma}}^{\mathbf{k}\alpha}( (E, \mathbf{t}) )
    &= \begin{pmatrix}
      \cos (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} & -\sin (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} \\
      \sin (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} & \cos (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} \\
    \end{pmatrix}
```


## References

```{bibliography}
:filter: docname in docnames
```
