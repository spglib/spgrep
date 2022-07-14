# On-the-fly irreps generation from solvable group chain

## Subduced and induced representations for solvable group

Consider projective irrep {math}`\Gamma` with factor system {math}`\mu` for solvable group {math}`G`,
```{math}
  R \phi_{j} = \sum_{i=1}^{d} \phi_{i} \Gamma(R)_{ij} \quad (R \in G, i = 1, \dots, d).
```
Let {math}`H` be invariant subgroup of {math}`G` such that {math}`G/H = \mathbb{Z}_{p}` ({math}`p` is prime and {math}`r^{p}=e`),
```{math}
  G = \coprod_{m=0}^{p-1} r^{m} H.
```

Let one of irrep decomposed from subduced representation {math}`\Gamma \downarrow H` be {math}`(\Delta, \mathrm{Span}\{ \psi_{i} \}_{i})`.
The following representations are also irrep over {math}`H`:
```{math}
  \Delta^{(m)}(S)
    &:= \frac{ \mu(S, r^{m}) }{ \mu(r^{m}, S_{m}) } \Delta(S_{m})
    \quad (S \in H, m = 0, \dots, p-1) \\
  S_{m}
    &:= r^{-m} S r^{m} (\in H)
```
There are two cases for {math}`\{ \Delta^{(m)} \}_{m=0, \dots, p-1}`:
1. They are mutually inequivalent: {math}`\Gamma \downarrow H \cong \sum_{m=0}^{p-1} \Delta^{(m)}`
2. They are all equivalent: {math}`\Gamma \downarrow H \cong p \Delta`

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

We sometimes need to restrict irrep under a vector space over {math}`\mathbb{R}` (instead of {math}`\mathbb{C}`), which is called physically irreducible representation [^SHW91].
When {math}`\Gamma` is real, **TODO**.

When {math}`\Gamma` is pseudo-real or not equivalent to {math}`\Gamma^{\ast}`, transform conjugated basis pair to real vectors by unitary matrix:
  ```{math}
    (\mathbf{v}_{i} \, \mathbf{v}_{i}^{\ast})
    \begin{pmatrix}
      \frac{1}{2} & -\frac{i}{2} \\
      \frac{1}{2} & \frac{i}{2}
    \end{pmatrix}
    = (\mathop{\mathrm{Re}} \mathbf{v}_{i} \, \mathop{\mathrm{Im}} \mathbf{v}_{i}).
  ```
  There is ambiguity to order vectors as {math}`(\mathbf{v}_{i} \, \mathbf{v}_{i}^{\ast})` or {math}`(\mathbf{v}_{i}^{\ast} \, \mathbf{v}_{i})`.
  We choose the former if {math}`\mathop{\mathrm{Im}} \mathbf{v}_{i}` is greater than {math}`\mathbf{0}` in lexicographic order and vice versa.

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

To symmetrize multi-dimensional irrep {math}`\Gamma^{(E)}`, consider its subduced representation over {math}`\{ g_{0}, g_{1}, g_{2} \} \cong C_{3}`.
```{math}
  \Gamma^{(E)} \downarrow C_{3} = \Delta^{(\,^{1}E)} + \Delta^{(\,^{2}E)}
```

|                       | {math}`g_{0}` | {math}`g_{1}` | {math}`g_{2}` |
| :-------------------- | ------------- | ------------- | ------------- |
| {math}`\Delta^{(\,^{1}E)}` | 1 | {math}`\omega=e^{i\frac{2}{3}\pi}` | {math}`\omega^{\ast}` |
| {math}`\Delta^{(\,^{2}E)}` | 1 | {math}`\omega^{\ast}` | {math}`\omega` |

When we choose {math}`g_{1}` as generator of {math}`C_{3}`, we select {math}`\Delta^{(E2)}` for generating induced representation [^footnote4].
Since {math}`\Delta^{(\,^{1}E)}` and {math}`\Delta^{(\,^{2}E)}` are inequivalent, induced representation {math}`\Gamma^{(E, sym)} = \Delta^{(\,^{2}E)} \uparrow C_{3v}` is irreducible.

[^footnote4]: Of course it is arbitrary which irreps to choose.
    Here, we choose irrep with the smallest angle of character of {math}`g_{1}` in {math}`(-\pi, \pi]`.
    Then, {math}`\mathrm{arg}(\chi^{(\,^{1}E)}(g_{1})) = \frac{2}{3}\pi` and {math}`\mathrm{arg}(\chi^{(\,^{2}E)}(g_{1})) = -\frac{2}{3}\pi`.
    Hence, we select {math}`\Delta^{(\,^{2}E)}`.

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

To symmetrize, consider subduced representation over {math}`\{ g_{0}, g_{2}, g_{4}, g_{6} \} \cong mm2` and {math}`\{ g_{0}, g_{4} \} \cong m`:
```{math}
  \Gamma^{(X_{1})} \downarrow mm2 &= \Delta^{(X_{1})} \\
  \Delta^{(X_{1})} \downarrow m &= \Theta^{(X_{1}, 1)} + \Theta^{(X_{1}, 2)}
```

|                             | {math}`g_{0}` | {math}`g_{4}` |
| :-------------------------- | ------------- | ------------: |
| {math}`\Theta^{(X_{1}, 1)}` | 1             | 1             |
| {math}`\Theta^{(X_{1}, 2)}` | 1             | -1            |

Since {math}`\Theta^{(X_{1}, 1)}` and {math}`\Theta^{(X_{1}, 2)}` are inequivalent, {math}`\Theta^{(X_{1}, 2)} \uparrow mm2` is irrep.


#### {math}`\mathcal{G} = Ia\overline{3}d` (No. 230) at {math}`H=(\frac{1}{2}\overline{\frac{1}{2}}\frac{1}{2})_{\mathrm{primitive}}`

(corresponding to {math}`G^{4}_{96}` in Ref. [^BC09])

```{math}
  \overline{\mathcal{G}}^{\mathbf{k}}| \cong m\overline{3}m,
  \quad
  \left| \overline{\mathcal{G}}^{\mathbf{k}} \right| = 48
```

```{math}
  \Gamma^{(\mathrm{reg})} = 2\Gamma^{(H_{1})} + 2\Gamma^{(H_{2})} + 6\Gamma^{(H_{3})}
```

[^SHW91]: Harold T. Stokes, Dorian M. Hatch, and James D. Wells, 
Phys. Rev. B 43, 11010 (1991).
