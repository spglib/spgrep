# On-the-fly irreps generation

## Regular representation

Given a finite group {math}`G = \{ R_{k} \}_{k=1}^{|G|}`, the regular representation of {math}`G` is defined as
```{math}
  \Gamma^{(\mathrm{reg})}(G_{k})_{ij}
    = \delta( G_{i}^{-1} G_{k} G_{j} ),
```
where {math}`\delta(\cdot)` takes one for identity operation and zero for others.

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
The factor system defined in "[](little_group.md)" satisfies this condition.

## Obtaining all irreps from (projective) regular representation

Inspired by Refs. [^Net73] and [^TVdV17].

For a finite group {math}`G` and its (projective) regular representation {math}`\Gamma^{(\mathrm{reg})}`,
```{math}
  \tilde{H}_{ij}
    = \sum_{R \in G} \sum_{k=1}^{|G|} \sum_{l=1}^{|G|} \Gamma^{(\mathrm{reg})}(R)_{ik} H_{kl} \Gamma^{(\mathrm{reg})}(R)_{jl}^{\ast}
```
is also Hermite and commute with {math}`\Gamma^{(\mathrm{reg})}`, where {math}`\mathbf{H}` is any Hermite matrix.
Consider the spectral decomposition of {math}`\tilde{\mathbf{H}}`,
```{math}
  \tilde{\mathbf{H}} = \sum_{\lambda} \lambda \sum_{n=1}^{ d_{\lambda} } \mathbf{v}_{\lambda n} \mathbf{v}_{\lambda n}^{\dagger} \\
  \lambda \in \mathbb{R}, \mathbf{v}_{\lambda n}^{\dagger} \mathbf{v}_{\lambda' n'} = \delta_{\lambda \lambda'} \delta_{n n'}.
```
Then, the regular representation is block-diagonalized with {math}`\mathbf{V} = ( \mathbf{v}_{\lambda_{1} 1} \dots \mathbf{v}_{\lambda_{1} d_{\lambda_{1}}} \dots )`,
```{math}
  \left[ \mathbf{V}^{\dagger} \mathbf{\Gamma}^{(\mathrm{reg})}(R) \mathbf{V} \right]_{\lambda n, \lambda' n'} (\lambda - \lambda') = 0
```

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

Note that Eqs. :eq:`all_irreps_dim_sum` also holds for projective representations.

## Working Example

### Crystallographic point group

Consider crystallographic point group {math}`3m` associated with {math}`P3m1` (No. 156).
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
\Gamma^{(\mathrm{ref})} = \Gamma^{(A_{1})} + \Gamma^{(A_{2})} + 2\Gamma^{(E)}
```

### Space group

Consider space group {math}`P4_{2}/mnm` (No. 136) at {math}`X=(0\frac{1}{2}0)`.

Consider irreps of space group {math}`\mathcal{G} = Ia\overline{3}d` (No. 230) at {math}`H=(\frac{1}{2}\overline{\frac{1}{2}}\frac{1}{2})_{\mathrm{primitive}}` (corresponding to {math}`G^{4}_{96}` in Ref. [^BC09]).

```{math}
  \overline{\mathcal{G}}^{\mathbf{k}}| \cong m\overline{3}m,
  \quad
  \left| \overline{\mathcal{G}}^{\mathbf{k}} \right| = 48
```

```{math}
  \Gamma^{(\mathrm{reg})} = 2\Gamma^{(H_{1})} + 2\Gamma^{(H_{2})} + 6\Gamma^{(H_{3})}
```


[^Net73]: N. Neto, Acta Cryst. A, 29(4) 464–472 (1973).
[^TVdV17]: John C. Thomas and Anton Van der Ven, J. Mech. Phys. Solids 107, 76–95, (2017).
[^BC09]: C. Bradley and A. P. Cracknell, The mathematical theory of symmetry in solids (Oxford, London, 2009).
