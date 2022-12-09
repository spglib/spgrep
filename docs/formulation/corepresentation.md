# Co-representation of magnetic point and space group

Refs: {cite}`doi:10.1126/sciadv.aat8685,wigner1959group`

## Anti-linear operator

Let $V$ be a vector space over $\mathbb{C}$ that operators in $M$ act.
An anti-unitary operator $a \in M$ satisfies
$$
    a (c_{1} \psi_{1} + c_{2} \psi_{2}) = c_{1}^{\ast} a\psi_{1} + c_{2}^{\ast} a\psi_{2}
$$
for $c_{1}, c_{2} \in \mathbb{C}$ and $\psi_{1}, \psi_{2} \in V$.


## Definition

Let $M$ be magnetic point group other than type I.
Let $D$ be a maximal point subgroup of $M$.
Then, $M$ can be decomposed as $M = D \sqcup a_{0} D$, where $a_{0}$ is an antisymmetry operation.

A co-representation $\Gamma$ gives linear or anti-linear operator for each element in $M$ while satisfying
$$
    \Gamma(u) \Gamma(u')        &= \omega(u, u') \Gamma(uu') \\
    \Gamma(u) \Gamma(a')        &= \omega(u, a') \Gamma(ua') \\
    \Gamma(a) \Gamma(u')^{\ast} &= \omega(a, u') \Gamma(au') \\
    \Gamma(a) \Gamma(a')^{\ast} &= \omega(a, a') \Gamma(aa'), \\
$$
where $u, u' \in D$ and $a, a' \in a_{0}D$.
The factor system $\omega$ satisfies a cocycle condition:
$$
    \omega(u, m) \omega(um, m') &= \omega(u, mm') \omega(m, m') \\
    \omega(a, m) \omega(am, m') &= \omega(a, mm') \omega(m, m')^{\ast} \\
$$
for $u \in D$, $a \in a_{0}D$, and $m, m \in M$ [^footnote1].

[^footnote1]: Consider $\Gamma(a)\Gamma(m)^{\ast}\Gamma(m')^{\ast}$ for $m \in D$ and $\Gamma(a)\Gamma(m)^{\ast}\Gamma(m')$ for $m \in a_{0}D$.

Co-representations $\Gamma$ and $\Gamma'$ are equivalent if an invertible matrix $\mathbf{T}$ exists such that
$$
    \Gamma'(u) &= \mathbf{T}^{-1} \Gamma(u) \mathbf{T} \quad (u \in D) \\
    \Gamma'(a) &= \mathbf{T}^{-1} \Gamma(a) \mathbf{T}^{\ast} \quad (a \in a_{0}D) \\
$$

## Frobenius-Schur indicator for co-representation

Let $(\Gamma, \mathrm{Span}_{\mathbb{C}} \{ \mathbf{\phi}_{i} \}_{i=1}^{d} )$ be one of **unitary** irreps for $D$ with factor system $\omega$.
Then, $\{ a_{0} \mathbf{v}_{i} \}_{i}$ also form irrep as
$$
    \overline{\Gamma}(u) &:= \frac{ \omega(u, a_{0}) }{ \omega( a_{0}, a_{0}^{-1} u a_{0} ) } \Gamma( a_{0}^{-1} u a_{0} )^{\ast} \quad (u \in D) \\
    u a_{0} \mathbf{\phi}_{j} &= \sum_{j} a_{0} \mathbf{\phi}_{j} \overline{\Gamma}(u)_{ij}.
$$

The following Frobenius-Schur indicator should be one of $\{ -1, 0, 1 \}$:
$$
    \xi^{\alpha} &:= \frac{1}{|D|} \sum_{ u \in D } \omega(a_{0}u, a_{0}u) \chi( (a_{0}u)^{2} ) \\
    \chi(u) &:= \mathrm{Tr} \Gamma(u) \quad (u \in D).
$$
This indicator can check if $\Gamma$ and $\overline{\Gamma}$ are equivalent.

### Case: $\xi^{\alpha} = 1$

In this case, $\Gamma$ and $\overline{\Gamma}$ are equivalent.
Let $\mathbf{U}$ be a **unitary** intertwiner between $\Gamma$ and $\overline{\Gamma}$:
$$
    \overline{\Gamma}(u) = \mathbf{U}^{-1} \Gamma(u) \mathbf{U} \quad (u \in D).
$$
Then, a transformed basis $\{ \mathbf{\psi}_{i} := \sum_{j} a_{0} \mathbf{\phi}_{j} [\mathbf{U}^{\dagger}]_{ji} \}_{i}$ also forms $\Gamma$.

A new basis $\{ \frac{1}{\sqrt{2}}(\phi_{i} + \psi_{i}) \}_{i}$ gives the following irrep [^derivation]:
$$
    \tilde{\mathbf{\Gamma}}(u) &= \mathbf{\Gamma}(u) \\
    \tilde{\mathbf{\Gamma}}(a_{0}) &= \mathbf{U} \\
    \tilde{\mathbf{\Gamma}}(a_{0}u) &= \omega(a_{0}, u)^{\ast} \tilde{\mathbf{\Gamma}}(a_{0}) \tilde{\mathbf{\Gamma}}(u) \\
    \mathbf{U}\mathbf{U}^{\ast} &= \omega(a_{0}, a_{0})\Gamma(a_{0}^{2}).
$$

[^derivation]: Be careful anti-linear operators,
$$
    a_{0} \psi_{i}
        &= \sum_{j} [U^{\dagger\ast}]_{ji} a_{0} (a_{0} \phi_{j}) \\
        &= \sum_{jk} \phi_{k} \omega(a_{0}, a_{0}) \Gamma(a_{0}^{2})_{kj} [U^{\dagger\ast}]_{ji} \\
        &= \sum_{k} \phi_{k} U_{ki}
$$

### Case: $\xi^{\alpha} = -1$

In this case, $\Gamma$ and $\overline{\Gamma}$ are equivalent.
Let $\mathbf{U}$ be a **unitary** intertwiner between $\Gamma$ and $\overline{\Gamma}$:
$$
    \overline{\Gamma}(u) &= \mathbf{U}^{-1} \Gamma(u) \mathbf{U} \quad (u \in D) \\
    \mathbf{U}\mathbf{U}^{\ast} &= -\omega(a_{0}, a_{0})\Gamma(a_{0}^{2}).
$$

We can take $(\mathbf{\phi}_{1}, \cdots, \mathbf{\phi}_{d}, \mathbf{\psi}_{1}, \cdots, \mathbf{\psi}_{d})$ as basis and they form irrep of $M$ as
$$
    \tilde{\mathbf{\Gamma}}(u) &=
        \begin{pmatrix}
            \mathbf{\Gamma}(u) & \mathbf{0} \\
            \mathbf{0} & \overline{\mathbf{\Gamma}}(u) \\
        \end{pmatrix} \\
    \tilde{\mathbf{\Gamma}}(a_{0}) &=
        \begin{pmatrix}
            \mathbf{0} & -\mathbf{U} \\
            \mathbf{U} & \mathbf{0} \\
        \end{pmatrix} \\
    \tilde{\mathbf{\Gamma}}(a_{0}u) &= \omega(a_{0}, u)^{\ast} \tilde{\mathbf{\Gamma}}(a_{0}) \tilde{\mathbf{\Gamma}}(u)
$$

### Case: $\xi^{\alpha} = 0$

In this case, $\Gamma$ and $\overline{\Gamma}$ are not equivalent.
We can take $(\mathbf{\phi}_{1}, \cdots, \mathbf{\phi}_{d}, a_{0}\mathbf{\phi}_{1}, \cdots, a_{0}\mathbf{\phi}_{d})$ as basis and they form irrep of $M$ as
$$
    \tilde{\mathbf{\Gamma}}(u) &=
        \begin{pmatrix}
            \mathbf{\Gamma}(u) & \mathbf{0} \\
            \mathbf{0} & \overline{\mathbf{\Gamma}}(u) \\
        \end{pmatrix} \\
    \tilde{\mathbf{\Gamma}}(a_{0}) &=
        \begin{pmatrix}
            \mathbf{0} & \omega(a_{0}, a_{0}) \mathbf{\Gamma}(a_{0}^{2}) \\
            \mathbf{1} & \mathbf{0} \\
        \end{pmatrix} \\
    \tilde{\mathbf{\Gamma}}(a_{0}u) &= \omega(a_{0}, u)^{\ast} \tilde{\mathbf{\Gamma}}(a_{0}) \tilde{\mathbf{\Gamma}}(u)
$$

(corep_spinor_factor_system)=
## Convention of anti-linear operators in spgrep

Let $M$ be magnetic point group other than type I.
Let $D$ be a maximal point subgroup of $M$.
Then, $M$ can be decomposed as $M = D \sqcup D a_{0}$, where $a_{0} = \mathbf{S}_{0} 1'$ is an antisymmetry operation.
We choose the following convention for choosing unitary or anti-unitary matrices for $M$:
$$
    D \ni \mathbf{S} &\mapsto \mathbf{U}(\mathbf{S}) \\
    D a_{0} \ni \mathbf{S}1' &\mapsto \mathbf{U}(\mathbf{S}) \left( -i \sigma_{y} \right) K,
    % a_{0} &\mapsto \mathbf{U}(\mathbf{S}_{0}) \left( -i \sigma_{y} \right) K \\
    % a_{0} D \ni a_{0} \mathbf{S} &\mapsto \mathbf{U}(\mathbf{S}_{0}) \left( -i \sigma_{y} \right) \mathbf{U}(\mathbf{S})^{\ast} K,
$$
where $\mathbf{U}: O(3) \to SU(2)$ is defined in {ref}`spinor_factor_system`.
$\sigma_{y}$ is the Pauli matrix.
$K$ is an anti-linear operator that takes complex conjugate.

## References

```{bibliography}
:filter: docname in docnames
```
