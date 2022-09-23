(intertwiner)=
# Numerically obtain intertwiner between equivalent representations

## Intertwiner between two projective representations of finite group

Let {math}`\Delta^{(0)}` and {math}`\Delta^{(1)}` be projective representations of a finite group {math}`H`.
An intertwiner of {math}`\Delta^{(0)}` and {math}`\Delta^{(1)}` is a matrix {math}`\mathbf{U}` satisfying
```{math}
    \mathbf{\Delta}^{(0)}(S) \mathbf{U} = \mathbf{U} \mathbf{\Delta}^{(1)}(S)
    \quad (\forall S \in H)
```

The following matrix is an intertwiner for these projective representations:
```{math}
  \mathbf{U} = \sum_{S \in H} \mathbf{\Delta}^{(0)}(S) \mathbf{B} \mathbf{\Delta}^{(1)}(S)^{-1},
```
where {math}`\mathbf{B}` is any matrix.

## Intertwiner between two space-group representations of space group {math}`\mathcal{G}`

Let {math}`\Gamma^{(\mathbf{k}, 0)}` and {math}`\Gamma^{(\mathbf{k}, 1)}` be representations of space group {math}`\mathcal{G}` with {math}`\mathbf{k}` vector.
Let {math}`\mathcal{T}` be a translational subgroup of {math}`\mathcal{G}`.
Consider coset representatives of little group of {math}`\mathbf{k}` over {math}`\mathcal{T}`:
```{math}
  \mathcal{G}^{\mathbf{k}} = \coprod_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } (\mathbf{S}_{i}, \mathbf{w}_{i}) \mathcal{T}.
```

Consider a similar matrix as the above finite-group case and take summation over lattice points beforehand:
```{math}
    \mathbf{U}
    &:= \frac{1}{N} \sum_{ \mathbf{t} } \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
            \mathbf{\Gamma}^{(\mathbf{k}, 0)}\left( (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) \right)
            \mathbf{B}
            \mathbf{\Gamma}^{(\mathbf{k}, 1)}\left( (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) \right)^{-1} \\
    &= \left(
            \frac{1}{N} \sum_{ \mathbf{t} } 1
       \right)
       \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
            \mathbf{\Gamma}^{(\mathbf{k}, 0)}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)
            \mathbf{B}
            \mathbf{\Gamma}^{(\mathbf{k}, 1)}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)^{-1} \\
    &= \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
            \mathbf{\Gamma}^{(\mathbf{k}, 0)}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)
            \mathbf{B}
            \mathbf{\Gamma}^{(\mathbf{k}, 1)}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)^{-1} \\
```
