# Irreps of space group

This page formulates little group, small representation, and irreducible representations (irreps) of space groups.

## Reciprocal lattice and transformation

Let {math}`\mathbf{A} = (\mathbf{a}_{1} \mathbf{a}_{2} \mathbf{a}_{3})` be basis of {math}`L_{\mathcal{T}}`, which is lattice of translational group {math}`{\mathcal{T}}`.
We choose basis of the (crystallographic) reciprocal vectors as
```{math}
  \mathbf{B} &= (\mathbf{b}_{1} \mathbf{b}_{2} \mathbf{b}_{3}) = A^{-\top} \\
  \mathbf{b}_{i} \cdot \mathbf{a}_{j} &= \delta_{ij}.
```

When basis vectors {math}`\mathbf{A}` is transformed to {math}`\mathbf{A}' := \mathbf{AP}`,
- its dual basis vectors in reciprocal space is transformed from {math}`\mathbf{B}=\mathbf{A}^{-\top}` to {math}`\mathbf{B}' := \mathbf{B}\mathbf{P}^{-\top}`.
- symmetry operation {math}`(\mathbf{R}, \mathbf{\tau})` is transformed to
    ```{math}
        (\mathbf{P}, \mathbf{0})^{-1} (\mathbf{R}, \mathbf{\tau}) (\mathbf{P}, \mathbf{0})
            = (\mathbf{P}^{-1}\mathbf{RP}, \mathbf{P}^{-1}\mathbf{\tau}).
    ```
- Coefficients of reciprocal vector {math}`\mathbf{k} = 2 \pi \mathbf{B} \mathbf{k}_{f}` is transformed to {math}`\mathbf{k}_{f}' := \mathbf{P}^{\top} \mathbf{k}_{f}`.

## Irreps of translation subgroup

Consider the following irreps of translation subgroup {math}`\mathcal{T}` of space group {math}`\mathcal{G}`,
```{math}
  \Gamma^{(\mathbf{k})} ((\mathbf{E}, \mathbf{t}))
  = \exp \left( -i \mathbf{k} \cdot \mathbf{t} \right)
  \quad ( \mathbf{t} \in L_{\mathcal{T}} ),
```
where {math}`L_{\mathcal{T}}` is a lattice formed by {math}`\mathcal{T}`.
Let {math}`L_{\mathcal{T}}^{\ast}` be a reciprocal lattice of {math}`L_{\mathcal{T}}`.
We can confine {math}`\mathbf{k}` within or the surface of a primitive cell of the reciprocal space because, for {math}`\mathbf{g} \in L_{\mathcal{T}}^{\ast}`, {math}`\Gamma^{(\mathbf{k})} = \Gamma^{(\mathbf{k} + \mathbf{g})}`.
Any Bloch function with {math}`\mathbf{k}`
```{math}
  \Psi_{\mathbf{k}}(\mathbf{r}) = e^{ i \mathbf{k} \cdot \mathbf{r} } u_{\mathbf{k}}(\mathbf{r}),
```
where {math}`u_{\mathbf{k}}(\mathbf{r})` is periodic with {math}`L_{\mathcal{T}}`, can be taken as (one-dimensional) basis of {math}`\Gamma^{\mathbf{k}}`:
```{math}
  (\mathbf{E}, \mathbf{t}) \Psi_{\mathbf{k}}(\mathbf{r})
    &= \Psi_{\mathbf{k}}( (\mathbf{E}, \mathbf{t})^{-1} \mathbf{r}) \\
    &= e^{-i\mathbf{k}\cdot\mathbf{t}} \Psi_{\mathbf{k}}(\mathbf{r}).
```

We write point group of space group {math}`\mathcal{G}` as {math}`\mathcal{P}`.
For {math}`(\mathbf{R}, \mathbf{v}) \in \mathcal{G}`, {math}`(\mathbf{R}, \mathbf{v}) \Psi_{\mathbf{k}}(\mathbf{r})` is Bloch function with {math}`\mathbf{Rk}`:
```{math}
:label: eq_point_group_in_bloch_function
  (\mathbf{E}, \mathbf{t}) (\mathbf{R}, \mathbf{v}) \Psi_{\mathbf{k}}(\mathbf{r})
  &= (\mathbf{R}, \mathbf{v}) (\mathbf{E}, \mathbf{R}^{-1}\mathbf{t}) \Psi_{\mathbf{k}}(\mathbf{r}) \nonumber \\
  &= \exp \left( -i \mathbf{k} \cdot \mathbf{R}^{-1} \mathbf{t} \right) (\mathbf{R}, \mathbf{v}) \Psi_{\mathbf{k}}(\mathbf{r}) \nonumber \\
  &= \exp \left( -i \mathbf{Rk} \cdot \mathbf{t} \right) (\mathbf{R}, \mathbf{v}) \Psi_{\mathbf{k}}(\mathbf{r}) \quad (\because \mathbf{R} \in O(3) ).
```
Be careful the basis for calculating inner product in Eq. {eq}`eq_point_group_in_bloch_function`, when you use crystallographic coordinates [^footnote1]:
```{math}
  \mathbf{v} &=: \mathbf{A} \mathbf{v}_{f} \\
  \mathbf{A}^{-1} \mathbf{R} \mathbf{A} &=: \mathbf{R}_{f} \quad \in \mathrm{SL}(3) \\
  (\mathbf{k}, \mathbf{R}^{-1} \mathbf{v})
    &= (\mathbf{R} \mathbf{k}, \mathbf{v}) \nonumber \\
    &= 2 \pi \left( ( \mathbf{A}^{\top} \mathbf{A}^{-1} ) \mathbf{R}_{f} ( \mathbf{A}^{\top} \mathbf{A}^{-1} )^{-1} \mathbf{k}_{f}, \mathbf{v}_{f} \right) \nonumber \\
    &= 2 \pi \left( \mathbf{R}_{f}^{\top} \mathbf{k}_{f}, \mathbf{v}_{f} \right)
    \quad (\because \mathbf{R}_{f}^{\top} = \mathbf{A}^{\top} \mathbf{R}^{-1} \mathbf{A}^{-\top} ).
```

[^footnote1]: We interchangeably denote an inner product of two vectors as {math}`\mathbf{a} \cdot \mathbf{b}` or {math}`(\mathbf{a}, \mathbf{b})`.

## Little group

Let {math}`\mathcal{P}` be point group of space group {math}`\mathcal{G}`.
We define little co-group of {math}`\mathbf{k}`,
```{math}
  \overline{\mathcal{G}}^{\mathbf{k}} = \left\{ \mathbf{R} \in \mathcal{P} \mid \mathbf{Rk} \equiv \mathbf{k} \right\},
```
where {math}`\equiv` is up to reciprocal vectors in {math}`L_{\mathcal{T}}^{\ast}`.
When we use crystallographic coordinates w.r.t. **primitive** basis vectors, the condition {math}`\mathbf{Rk} \equiv \mathbf{k}` is written as
```{math}
    \mathbf{Rk} \equiv \mathbf{k}
    \Longleftrightarrow
    \mathbf{R}_{f}^{\top} \mathbf{k}_{f} - \mathbf{k}_{f} \in \mathbb{Z}^{3}.
```

Because lattice {math}`L_{\mathcal{T}}` is invariant with {math}`\mathbf{R} \in \mathcal{P}`,
```{math}
  \mathbf{R} \mathbf{k}_{1} \equiv \mathbf{k}_{2}
  \Rightarrow
  \overline{\mathcal{G}}^{\mathbf{k}_{2}} = \mathbf{R} \overline{\mathcal{G}}^{\mathbf{k}_{1}} \mathbf{R}^{-1}.
```
When we decompose point group {math}`\mathcal{P}` as
```{math}
  \mathcal{P} = \coprod_{i=1}^{q} \mathbf{P}_{i} \overline{\mathcal{G}}^{\mathbf{k}_{1}},
```
{math}`\mathbf{P}_{i}\mathbf{k}_{1} \equiv \mathbf{k}_{i}` with star of {math}`\mathbf{k}_{1}`, {math}`\left\{ \mathbf{k}_{1}, \dots, \mathbf{k}_{q} \right\}`.

We define little group (of the first kind) of {math}`\mathbf{k}` as
```{math}
  \mathcal{G}^{\mathbf{k}} = \coprod_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } (\mathbf{S}_{i}, \mathbf{w}_{i}) \mathcal{T},
```
where {math}`\mathcal{G} = \coprod_{i} (\mathbf{S}_{i}, \mathbf{w}_{i}) \mathcal{T}`.
We can decompose space group {math}`\mathcal{G}` as
```{math}
  \mathcal{G} = \coprod_{i=1}^{q} (\mathbf{P}_{i}, \mathbf{x}_{i}) \mathcal{G}^{\mathbf{k}_{1}}.
```
Then, for Bloch function {math}`\Psi_{\mathbf{k}_{1}}(\mathbf{r})`, {math}`(\mathbf{P}_{i}, \mathbf{x}_{i}) \mathcal{G}^{\mathbf{k}_{1}} \Psi_{\mathbf{k}_{1}}(\mathbf{r})` is Bloch function for {math}`\mathbf{k}_{i}`.
Thus, our task to calculate irreps of {math}`\mathcal{G}` is boiled down to calculating irreps of little group {math}`\mathcal{G}^{\mathbf{k}}`, which is called small representation.

## Projective representation

A projective representation of group {math}`G` is a non-singular matrix function {math}`\Delta` on {math}`G` if it satisfied the following conditions:
For each group product {math}`G_{k} = G_{i}G_{j}`, there exists scalar function {math}`\mu` (factor system) such that
```{math}
  \Delta(G_{i}) \Delta(G_{j}) = \mu(G_{i}, G_{j}) \Delta(G_{k}),
```
and
```{math}
  \mu(G_{i}, G_{j}G_{k}) \mu(G_{j}, G_{k}) = \mu(G_{i}G_{j}, G_{k}) \mu(G_{i}, G_{j})
  \quad \mbox{(cocycle condition)}.
```

## Small representation

Let {math}`\Gamma^{\mathbf{k}}_{p}` be irrep of {math}`\mathcal{G}^{\mathbf{k}}`, for {math}`\mathbf{S}_{i} \mathbf{S}_{j} = \mathbf{S}_{k}`,
```{math}
  \Gamma^{\mathbf{k}}_{p}((\mathbf{S}_{i}, \mathbf{w}_{i})) \Gamma^{\mathbf{k}}_{p}((\mathbf{S}_{j}, \mathbf{w}_{j}))
    = \exp \left( -i \mathbf{k} \cdot ( \mathbf{w}_{i} + \mathbf{S}_{i} \mathbf{w}_{j} - \mathbf{w}_{k} ) \right) \Gamma^{\mathbf{k}}_{p}((\mathbf{S}_{k}, \mathbf{w}_{k})).
```
Then, Simplifying by
```{math}
  \Gamma^{\mathbf{k}}_{p}((\mathbf{R}, \mathbf{v})) &=: \exp \left( -i \mathbf{k} \cdot \mathbf{v} \right) D^{\mathbf{k}}_{p}((\mathbf{R}, \mathbf{v})) \\
  D^{\mathbf{k}}_{p}((\mathbf{S}_{i}, \mathbf{w}_{i})) D^{\mathbf{k}}_{p}((\mathbf{S}_{j}, \mathbf{w}_{j}))
    &= \exp \left( -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} \right) D^{\mathbf{k}}_{p}((\mathbf{S}_{k}, \mathbf{w}_{k})),
```
where
```{math}
  \mathbf{g}_{i}
    &= \mathbf{S}_{i}^{-1} \mathbf{k} - \mathbf{k} \\
    &= 2 \pi \mathbf{B} \left( \mathbf{S}_{i, f}^{\top} \mathbf{k}_{f} - \mathbf{k}_{f} \right)
        \quad (\mathbf{S}_{i, f} = \mathbf{A} \mathbf{S}_{i} \mathbf{A}^{-1}).
```

**If we choose origin such that {math}`(\mathbf{E}, \mathbf{t}) \in \mathcal{G}`**, {math}`D^{\mathbf{k}}_{p}((\mathbf{E}, \mathbf{t})) = \mathbf{1}` and {math}`\exp \left( -i \mathbf{g}_{i} \cdot (\mathbf{w}_{j} + \mathbf{t}) \right) =\exp \left( -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} \right) ` hold.
Then, it is sufficient to consider on {math}`\overline{\mathcal{G}}^{\mathbf{k}}`,
```{math}
  D^{\mathbf{k}}_{p}(\mathbf{S}_{i}) D^{\mathbf{k}}_{p}(\mathbf{S}_{j})
    = \exp \left( -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} \right) D^{\mathbf{k}}_{p}(\mathbf{S}_{k})
```
Here, {math}`\mu_{\mathbf{k}}(S_{i}, S_{j}) := \exp \left( -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} \right)` holds the cocycle condition and {math}`\mu_{\mathbf{k}}(\mathbf{E}, \mathbf{E}) = 1`.
Now, our task is to enumerate projective irreps of little co-group {math}`\overline{\mathcal{G}}^{\mathbf{k}}` with factor system {math}`\mu_{\mathbf{k}}`.

When we use crystallographic coordinates, the factor system is written as
```{math}
    \mu_{\mathbf{k}}(\mathbf{S}_{i}, \mathbf{S}_{j})
        &= \exp \left( -i (\mathbf{g}_{i}, \mathbf{w}_{j}) \right) \\
        &= \exp \left( -2 \pi i (\mathbf{S}_{i, f}^{\top} \mathbf{k}_{f} - \mathbf{k}_{f}, \mathbf{w}_{j, f}) \right)
            \quad (\mathbf{w}_{j} =: \mathbf{A} \mathbf{w}_{j, f}).
```
