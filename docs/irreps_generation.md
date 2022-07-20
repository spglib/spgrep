# On-the-fly irreps generation from solvable group chain

## Subduced and induced representations for solvable group

Let {math}`H` be invariant subgroup of {math}`G` such that {math}`G/H = \mathbb{Z}_{p}` ({math}`p` is prime and {math}`r^{p}=e`),
```{math}
  G = \coprod_{m=0}^{p-1} r^{m} H.
```
Consider projective irrep {math}`\Delta` with factor system {math}`\mu` for {math}`H`,
```{math}
  S \psi_{j} = \sum_{i=1}^{d} \psi_{i} \Delta(S)_{ij} \quad (R \in H, i = 1, \dots, d).
```
The following representations are also irrep over {math}`H`:
```{math}
  \Delta^{(m)}(S)
    &:= \frac{ \mu(S, r^{m}) }{ \mu(r^{m}, S_{m}) } \Delta(S_{m})
    \quad (S \in H, m = 0, \dots, p-1) \\
  S_{m}
    &:= r^{-m} S r^{m} (\in H)
```
There are two cases for {math}`\{ \Delta^{(m)} \}_{m=0, \dots, p-1}`:
1. They are mutually inequivalent: {math}`\Delta \uparrow G \downarrow H \cong \sum_{m=0}^{p-1} \Delta^{(m)}`
2. They are all equivalent: {math}`\Delta \uparrow G \downarrow H \cong p \Delta`

### Case-1: conjugated irreps are mutually inequivalent

For case-1, the induced representation {math}`\Delta \uparrow G` is irrep [^footnote1] [^footnote2]:
```{math}
  \Delta \uparrow G (r^{m} S)_{i:, j:}
    &= \mathbb{I}[i \equiv j + m] \frac{ \mu(r^{m}S, r^{j}) }{ \mu(r^{i}, S_{j} ) } \Delta( S_{j} )
    \quad (S \in H).\\
  \Delta \uparrow G (S)
    &= \begin{pmatrix}
      \Delta^{(0)}(S) & & & \\
      & \Delta^{(1)}(S) & & \\
      & & \ddots & \\
      & & & \Delta^{(p-1)}(S)
    \end{pmatrix}
    \quad (S \in H) \\
  \Delta \uparrow G (r)
    &= \begin{pmatrix}
      \mathbf{0} & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} & \frac{\mu(r, r^{p-1})}{\mu(E,E)} \mathbf{1} \\
      \frac{\mu(r, E)}{\mu(r,E)} \mathbf{1} & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} \\
      \mathbf{0} & \frac{\mu(r, r)}{\mu(r^{2},E)} \mathbf{1} & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} \\
      \vdots     & \vdots     & \vdots     &        & \vdots     & \vdots \\
      \mathbf{0} & \mathbf{0} & \mathbf{0} & \cdots & \frac{\mu(r, r^{p-2})}{\mu(r^{p-1},E)} \mathbf{1} & \mathbf{0}
    \end{pmatrix}.
```

### Case-2: conjugated irreps are equivalent

[^footnote1]: This follows from Frobenius reciprocity theorem. In this case, {math}`\Delta` appears only once in {math}`\Gamma \downarrow H`.
Thus, {math}`\Gamma` appears only once in {math}`\Delta \uparrow G`.

[^footnote2]: We can choose factor system to satisfy {math}`\mu(E, R) = \mu(R, E) = 1 (\forall R \in G)`, and then {math}`\Delta^{(0)} = \Delta`

For case-2, the induced representation {math}`\Delta \uparrow G` is reducible.
Let one of intertwiner between {math}`\Delta^{(0)}` and {math}`\Delta^{(1)}` be {math}`\mathbf{U}`,
```{math}
  \mathbf{\Delta}^{(0)}(S) \mathbf{U} = \mathbf{U} \mathbf{\Delta}^{(1)}(S)
  \quad (\forall S \in H).
```
The intertwiner is unique up to scalar multiplication: if {math}`\mathbf{U}` and {math}`\mathbf{U}'` are intertwiner between {math}`\Delta^{(0)}` and {math}`\Delta^{(1)}`, {math}`\mathbf{U}^{-1}\mathbf{U}'` should be written as {math}`c\mathbf{I}` with some complex number {math}`c` from Schur's lemma.
The following matrix is an intertwiner for projective representations, {math}`\Delta^{(0)}` and {math}`\Delta^{(1)}`,
```{math}
  \mathbf{U} = \sum_{S \in H} \mathbf{\Delta}^{(0)}(S) \mathbf{B} \mathbf{\Delta}^{(1)}(S)^{-1},
```
where {math}`\mathbf{B}` is any matrix.
We scale {math}`\mathbf{U}` such that {math}`\mathbf{U}^{p} = \mathbf{1}`.

The induced representation {math}`\Delta \uparrow G` is decomposed to {math}`p` irreps {math}`\{ \Delta_{q} \}_{q=0}^{p-1}`,
```{math}
  \mathbf{\Delta_{q}}(S) &= \mathbf{\Delta}(S) \quad (S \in H) \\
  \mathbf{\Delta_{q}}(r) &= \frac{1}{\omega_{q}} \mathbf{U} \\
  \omega_{q}
    &:= \left( \frac{\prod_{m=1}^{p-1} \mu(r, r^{m})}{\mu(E, E)} \right)^{\frac{1}{p}} \exp \left( \frac{2 \pi i q}{p} \right).
```

## Reality of (projective) irrep

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

### (1) {math}`\Gamma` is real

{cite}`Stokes:pc5025,Inui1996-et`

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

## Working example

### Crystallographic point group

#### {math}`C_{3v}` associated with {math}`P3m1` (No. 156)

The matrix representation of this crystallographic point group is

```{math}
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
```

```{math}
\Gamma^{(\mathrm{reg})} = \Gamma^{(A_{1})} + \Gamma^{(A_{2})} + 2\Gamma^{(E)}
```

### Space group

#### {math}`P4_{2}/mnm` (No. 136) at {math}`X=(0\frac{1}{2}0)`

The little co-group is
```{math}
  \overline{\mathcal{G}}^{\mathbf{k}}
  &= \left\{
    g_{0} = e,
    g_{1},
    g_{2},
    g_{3} = g_{1} g_{2},
    g_{4},
    g_{5} = g_{1} g_{4},
    g_{6} = g_{2} g_{4},
    g_{7} = g_{1} g_{2} g_{4}
  \right\}
  \cong mmm, \\
  \mathrm{where} \quad
  g_{1} &= \mathrm{diag}(-1, -1, -1) \\
  g_{2} &= \mathrm{diag}(-1, -1, 1) \\
  g_{4} &= \mathrm{diag}(1, -1, -1).
```

```{math}
  \Gamma^{(\mathrm{reg})} = 2 \Gamma^{(X_{1})} + 2 \Gamma^{(X_{2})}
```

#### {math}`\mathcal{G} = Ia\overline{3}d` (No. 230) at {math}`H=(\frac{1}{2}\overline{\frac{1}{2}}\frac{1}{2})_{\mathrm{primitive}}`

(corresponding to {math}`G^{4}_{96}` in Ref. {cite}`Bradley2009-ze`)

```{math}
  \overline{\mathcal{G}}^{\mathbf{k}}| \cong m\overline{3}m,
  \quad
  \left| \overline{\mathcal{G}}^{\mathbf{k}} \right| = 48
```

```{math}
  \Gamma^{(\mathrm{reg})} = 2\Gamma^{(H_{1})} + 2\Gamma^{(H_{2})} + 6\Gamma^{(H_{3})}
```

## References

```{bibliography}
:filter: docname in docnames
```
