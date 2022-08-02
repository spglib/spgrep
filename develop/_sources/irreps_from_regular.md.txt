# On-the-fly irreps generation from regular representation

This page describes an algorithm for generating all irreps by random matrix.

## Regular representation

Given a finite group {math}`G = \{ R_{k} \}_{k=1}^{|G|}`, the regular representation of {math}`G` is defined as
```{math}
  \Gamma^{(\mathrm{reg})}(G_{k})_{ij}
    = \delta( G_{i}, G_{k} G_{j} ),
```
where {math}`\delta(\cdot, \cdot)` takes one iff two elements are equal.

For any factor system {math}`\mu`, the following representation holds the orthogonality theorem,
```{math}
  \Delta^{(\mathrm{reg})}(G_{k})_{ij}
    = \delta(G_{i}, G_{k}G_{j}) \mu(G_{k}, G_{j}).
```
We refer this representation to projective regular representation.

The projective regular representation is a unitary representation when absolute of its factor system is one:
```{math}
  \sum_{l} \Delta^{(\mathrm{reg})}(G_{k})_{il} \Delta^{(\mathrm{reg})}(G_{k})_{jl}^{\ast}
  = \delta_{ij} \left| \mu(G_{k}, G_{k}^{-1}G_{i}) \right|^{2}.
```
The factor system defined in "[](spacegroup_irreps.md)" satisfies this condition.

## Obtaining all irreps from (projective) regular representation

Refs. {cite}`THOMAS201776,doi:10.1137/090779966,Po2017`

For a finite group {math}`G` and its (projective) regular representation {math}`\Gamma^{(\mathrm{reg})}`,
```{math}
  \tilde{\mathbf{H}}
    = \sum_{R \in G} \mathbf{\Gamma}^{(\mathrm{reg})}(R) \mathbf{H} \mathbf{\Gamma}^{(\mathrm{reg})}(R^{-1})
```
is also Hermite and commute with {math}`\Gamma^{(\mathrm{reg})}` [^footnote0], where {math}`\mathbf{H}` is any Hermite matrix.
Consider the spectral decomposition of {math}`\tilde{\mathbf{H}}`,
```{math}
  \tilde{\mathbf{H}} = \sum_{\lambda} \lambda \sum_{n=1}^{ d_{\lambda} } \mathbf{v}_{\lambda n} \mathbf{v}_{\lambda n}^{\dagger} \\
  \lambda \in \mathbb{R}, \mathbf{v}_{\lambda n}^{\dagger} \mathbf{v}_{\lambda' n'} = \delta_{\lambda \lambda'} \delta_{n n'}.
```
Then, the regular representation is block-diagonalized with {math}`\mathbf{V} = ( \mathbf{v}_{\lambda_{1} 1} \dots \mathbf{v}_{\lambda_{1} d_{\lambda_{1}}} \dots )`,
```{math}
  \left[ \mathbf{V}^{\dagger} \mathbf{\Gamma}^{(\mathrm{reg})}(R) \mathbf{V} \right]_{\lambda n, \lambda' n'} (\lambda - \lambda') = 0
```

[^footnote0]: The derivation is as follows:
    ```{math}
      \mathbf{\Delta}^{(\mathrm{reg})}(S) \mathbf{U} \mathbf{\Delta}^{(\mathrm{reg})}(S)^{-1}
      &= \mathbf{\Delta}^{(\mathrm{reg})}(S) \sum_{P \in H} \mathbf{\Delta}^{(\mathrm{reg})}(P) \mathbf{\Delta}^{(\mathrm{reg})}(P)^{-1} \mathbf{\Delta}^{(\mathrm{reg})}(S)^{-1} \\
      &= \sum_{P \in H} \mu(S, P) \mathbf{\Delta}^{(\mathrm{reg})}(SP) \left( \mu(S, P) \mathbf{\Delta}^{(\mathrm{reg})}(SP) \right)^{-1} \\
      &= \sum_{P \in H} \mathbf{\Delta}^{(\mathrm{reg})}(P) \mathbf{\Delta}^{(\mathrm{reg})}(P)^{-1} \\
      &= \mathbf{U}
    ```
    Be careful that {math}`\mathbf{\Delta}^{(\mathrm{reg})}(S^{-1}) = \mu(S, S^{-1}) \mu(E, E) \mathbf{\Delta}^{(\mathrm{reg})}(S)^{-1}`.


{math}`\mathbf{V}_{\lambda} = ( \mathbf{v}_{\lambda 1} \dots \mathbf{v}_{\lambda_{1} d_{\lambda}} ) \in \mathbb{C}^{ |G| \times d_{\lambda} }` forms Irrep,
```{math}
  \mathbf{\Gamma}^{(\lambda)}(R)
    = \mathbf{V}_{\lambda}^{\dagger} \mathbf{\Gamma}^{(\mathrm{reg})}(R) \mathbf{V}_{\lambda}
    \quad
    \in \mathbb{C}^{ d_{\lambda} \times d_{\lambda} }
```
The regular representation contains all Irreps of the group up to unitary transformation as
```{math}
  \Gamma^{(\mathrm{reg})} = \sum_{\alpha} d_{\alpha} \Gamma^{(\alpha)},
```
where {math}`d_{\alpha}` is dimension of Irrep labeled as {math}`\alpha`, and the summation is taken over all Irreps up to unitary transformation.
Thus, this procedure produces all Irreps exhaustively.
For an infinite group, the above procedure also holds replacing regular representation into projective regular representation.

Also, we can check the obtained Irreps are enough by checking the following equality

```{math}
---
label: all_irreps_dim_sum
---
    \sum_{\alpha} d_{\alpha}^{2} &= |G| \\
    \sum_{ g \in G } | \chi^{(\alpha)}(g) |^2 &= |G| \\
    \sum_{g \in G} \chi^{(\alpha)}(g)^{\ast} \chi^{(\beta)}(g) &= |G| \delta_{\alpha, \beta},
```
where {math}`\chi^{(\alpha)}` is character of irrep {math}`\Gamma^{(\alpha)}`.

Note that Eqs. {eq}`all_irreps_dim_sum` also holds for projective representations.

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
