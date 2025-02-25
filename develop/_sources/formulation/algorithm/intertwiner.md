(intertwiner)=
# Numerically obtain intertwiner between equivalent representations

## Intertwiner between two projective representations of finite group

Let $\Delta^{(0)}$ and $\Delta^{(1)}$ be projective representations of a finite group $H$.
An intertwiner of $\Delta^{(0)}$ and $\Delta^{(1)}$ is a matrix $\mathbf{U}$ satisfying
$$
    \mathbf{\Delta}^{(0)}(S) \mathbf{U} = \mathbf{U} \mathbf{\Delta}^{(1)}(S)
    \quad (\forall S \in H)
$$
The intertwiner is unique up to scalar multiplication: if $\mathbf{U}$ and $\mathbf{U}'$ are intertwiner between $\Delta^{(0)}$ and $\Delta^{(1)}$, $\mathbf{U}^{-1}\mathbf{U}'$ should be written as $c\mathbf{I}$ with some complex number $c$ from Schur's lemma.
In particular, when $\Delta^{(0)}$ and $\Delta^{(1)}$ are unitary irreps, we can choose $\mathbf{U}$ as a unitary matrix [^footnote_unitary].

[^footnote_unitary]: [Equivalent unitary representations are unitarily equivalent](http://www.nathankarst.com/blog/equivalent-unitary-representations-are-unitarily-equivalent)

The following matrix is an intertwiner for these projective representations:
$$
  \mathbf{U} = \sum_{S \in H} \mathbf{\Delta}^{(0)}(S) \mathbf{B} \mathbf{\Delta}^{(1)}(S)^{-1},
$$
where $\mathbf{B}$ is any matrix.

## Intertwiner between two space-group representations of space group $\mathcal{G}$

Let $\Gamma^{(\mathbf{k}, 0)}$ and $\Gamma^{(\mathbf{k}, 1)}$ be representations of space group $\mathcal{G}$ with $\mathbf{k}$ vector.
Let $\mathcal{T}$ be a translational subgroup of $\mathcal{G}$.
Consider coset representatives of little group of $\mathbf{k}$ over $\mathcal{T}$:
$$
  \mathcal{G}^{\mathbf{k}} = \coprod_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } (\mathbf{S}_{i}, \mathbf{w}_{i}) \mathcal{T}.
$$

Consider a similar matrix as the above finite-group case and take summation over lattice points beforehand:
$$
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
$$
