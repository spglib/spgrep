# On-the-fly irreps generation from regular representation

Understand how diagonalizing a matrix that commutes with the regular representation produces irreps.

!!! info "Before you start"

    Read the [algorithm overview](algorithm.md) and [projective representations](../irreps/spacegroup_irreps.md#projective-representation).
    The worked cases below illustrate the mathematics. For runnable Python, use the [quickstart](../../quickstart.md).

This page describes an algorithm for generating all irreps by random matrix.

## Regular representation

Given a finite group $G = \{ R_{k} \}_{k=1}^{|G|}$, the regular representation of $G$ is defined as

$$
\begin{aligned}
\Gamma^{(\mathrm{reg})}(G_{k})_{ij}
  = \delta( G_{i}, G_{k} G_{j} ),
\end{aligned}
$$

where $\delta(\cdot, \cdot)$ takes one iff two elements are equal.

For any factor system $\mu$, the following representation holds the orthogonality theorem,

$$
\begin{aligned}
\Delta^{(\mathrm{reg})}(G_{k})_{ij}
  = \delta(G_{i}, G_{k}G_{j}) \mu(G_{k}, G_{j}).
\end{aligned}
$$

We refer this representation to projective regular representation.

The projective regular representation is a unitary representation when absolute of its factor system is one:

$$
\begin{aligned}
\sum_{l} \Delta^{(\mathrm{reg})}(G_{k})_{il} \Delta^{(\mathrm{reg})}(G_{k})_{jl}^{\ast}
= \delta_{ij} \left| \mu(G_{k}, G_{k}^{-1}G_{i}) \right|^{2}.
\end{aligned}
$$

The factor system defined in [Irreps of space group](../irreps/spacegroup_irreps.md#space_group_irreps) satisfies this condition.

## Obtaining all irreps from (projective) regular representation

Refs. [Thomas and {Van der Ven}, 2017](../../references.md#THOMAS201776), [Maehara and Murota, 2011](../../references.md#doi:10.1137/090779966), [Po et al., 2017](../../references.md#Po2017)

For a finite group $G$ and its (projective) regular representation $\Gamma^{(\mathrm{reg})}$,

$$
\begin{aligned}
\tilde{\mathbf{H}}
  = \sum_{R \in G} \mathbf{\Gamma}^{(\mathrm{reg})}(R) \mathbf{H} \mathbf{\Gamma}^{(\mathrm{reg})}(R^{-1})
\end{aligned}
$$

is also Hermite and commute with $\Gamma^{(\mathrm{reg})}$ [^footnote0], where $\mathbf{H}$ is any Hermite matrix.
Consider the spectral decomposition of $\tilde{\mathbf{H}}$,

$$
\begin{aligned}
\tilde{\mathbf{H}} = \sum_{\lambda} \lambda \sum_{n=1}^{ d_{\lambda} } \mathbf{v}_{\lambda n} \mathbf{v}_{\lambda n}^{\dagger} \\
\lambda \in \mathbb{R}, \mathbf{v}_{\lambda n}^{\dagger} \mathbf{v}_{\lambda' n'} = \delta_{\lambda \lambda'} \delta_{n n'}.
\end{aligned}
$$

Then, the regular representation is block-diagonalized with $\mathbf{V} = ( \mathbf{v}_{\lambda_{1} 1} \dots \mathbf{v}_{\lambda_{1} d_{\lambda_{1}}} \dots )$,

$$
\begin{aligned}
\left[ \mathbf{V}^{\dagger} \mathbf{\Gamma}^{(\mathrm{reg})}(R) \mathbf{V} \right]_{\lambda n, \lambda' n'} (\lambda - \lambda') = 0
\end{aligned}
$$

[^footnote0]: The derivation is as follows:

    $$
    \begin{aligned}
    \mathbf{\Delta}^{(\mathrm{reg})}(S) \mathbf{U} \mathbf{\Delta}^{(\mathrm{reg})}(S)^{-1}
    &= \mathbf{\Delta}^{(\mathrm{reg})}(S) \sum_{P \in H} \mathbf{\Delta}^{(\mathrm{reg})}(P) \mathbf{\Delta}^{(\mathrm{reg})}(P)^{-1} \mathbf{\Delta}^{(\mathrm{reg})}(S)^{-1} \\
    &= \sum_{P \in H} \mu(S, P) \mathbf{\Delta}^{(\mathrm{reg})}(SP) \left( \mu(S, P) \mathbf{\Delta}^{(\mathrm{reg})}(SP) \right)^{-1} \\
    &= \sum_{P \in H} \mathbf{\Delta}^{(\mathrm{reg})}(P) \mathbf{\Delta}^{(\mathrm{reg})}(P)^{-1} \\
    &= \mathbf{U}
    \end{aligned}
    $$

    Be careful that $\mathbf{\Delta}^{(\mathrm{reg})}(S^{-1}) = \mu(S, S^{-1}) \mu(E, E) \mathbf{\Delta}^{(\mathrm{reg})}(S)^{-1}$.

$\mathbf{V}_{\lambda} = ( \mathbf{v}_{\lambda 1} \dots \mathbf{v}_{\lambda_{1} d_{\lambda}} ) \in \mathbb{C}^{ |G| \times d_{\lambda} }$ forms Irrep,

$$
\begin{aligned}
\mathbf{\Gamma}^{(\lambda)}(R)
  = \mathbf{V}_{\lambda}^{\dagger} \mathbf{\Gamma}^{(\mathrm{reg})}(R) \mathbf{V}_{\lambda}
  \quad
  \in \mathbb{C}^{ d_{\lambda} \times d_{\lambda} }
\end{aligned}
$$

The regular representation contains all Irreps of the group up to unitary transformation as

$$
\begin{aligned}
\Gamma^{(\mathrm{reg})} = \sum_{\alpha} d_{\alpha} \Gamma^{(\alpha)},
\end{aligned}
$$

where $d_{\alpha}$ is dimension of Irrep labeled as $\alpha$, and the summation is taken over all Irreps up to unitary transformation.
Thus, this procedure produces all Irreps exhaustively.
For an infinite group, the above procedure also holds replacing regular representation into projective regular representation.

Also, we can check the obtained Irreps are enough by checking the following equality

<a id="eqn:all_irreps_dim_sum"></a>

$$
\begin{aligned}
\sum_{\alpha} d_{\alpha}^{2} &= |G| \\
\sum_{ g \in G } | \chi^{(\alpha)}(g) |^2 &= |G| \\
\sum_{g \in G} \chi^{(\alpha)}(g)^{\ast} \chi^{(\beta)}(g) &= |G| \delta_{\alpha, \beta},
\end{aligned}
\tag{1}
$$

where $\chi^{(\alpha)}$ is character of irrep $\Gamma^{(\alpha)}$.
Note that Eqs. [(1)](#eqn:all_irreps_dim_sum) also holds for projective representations.

## Working example

### Crystallographic point group

#### $C_{3v}$ associated with $P3m1$ (No. 156)

The matrix representation of this crystallographic point group is

$$
\begin{aligned}
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
\end{aligned}
$$

$$
\begin{aligned}
\Gamma^{(\mathrm{reg})} = \Gamma^{(A_{1})} + \Gamma^{(A_{2})} + 2\Gamma^{(E)}
\end{aligned}
$$

### Space group

#### $P4_{2}/mnm$ (No. 136) at $X=(0\frac{1}{2}0)$

The little co-group is

$$
\begin{aligned}
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
\end{aligned}
$$

$$
\begin{aligned}
\Gamma^{(\mathrm{reg})} = 2 \Gamma^{(X_{1})} + 2 \Gamma^{(X_{2})}
\end{aligned}
$$

#### $\mathcal{G} = Ia\overline{3}d$ (No. 230) at $H=(\frac{1}{2}\overline{\frac{1}{2}}\frac{1}{2})_{\mathrm{primitive}}$

(corresponding to $G^{4}_{96}$ in Ref. [Bradley and Cracknell, 2009](../../references.md#Bradley2009-ze))

$$
\begin{aligned}
\overline{\mathcal{G}}^{\mathbf{k}}| \cong m\overline{3}m,
\quad
\left| \overline{\mathcal{G}}^{\mathbf{k}} \right| = 48
\end{aligned}
$$

$$
\begin{aligned}
\Gamma^{(\mathrm{reg})} = 2\Gamma^{(H_{1})} + 2\Gamma^{(H_{2})} + 6\Gamma^{(H_{3})}
\end{aligned}
$$

## Next steps

Compare the [subgroup-chain construction](irreps_from_chain.md), or consult the [enumeration reference](../../api/api_rep.md#spgrep.rep.enumerate) for the implementation.
