# On-the-fly irreps generation from regular representation

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

Refs. {cite}`Neto:a09740,THOMAS201776,doi:10.1137/090779966,Po2017`

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

## References

```{bibliography}
:filter: docname in docnames
```
