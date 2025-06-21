# On-the-fly irreps generation from regular representation

This page describes an algorithm for generating all irreps by random matrix.

## Regular representation

Given a finite group $G = { R\_{k} }_{k=1}^{|G|}$, the regular representation of $G$ is defined as
$$
\\Gamma^{(\\mathrm{reg})}(G_{k})_{ij}
= \\delta( G_{i}, G\_{k} G\_{j} ),
$$
where $\\delta(\\cdot, \\cdot)$ takes one iff two elements are equal.

For any factor system $\\mu$, the following representation holds the orthogonality theorem,
$$
\\Delta^{(\\mathrm{reg})}(G\_{k})_{ij}
= \\delta(G_{i}, G\_{k}G\_{j}) \\mu(G\_{k}, G\_{j}).
$$
We refer this representation to projective regular representation.

The projective regular representation is a unitary representation when absolute of its factor system is one:
$$
\\sum\_{l} \\Delta^{(\\mathrm{reg})}(G\_{k})_{il} \\Delta^{(\\mathrm{reg})}(G_{k})_{jl}^{\\ast}
= \\delta_{ij} \\left| \\mu(G\_{k}, G\_{k}^{-1}G\_{i}) \\right|^{2}.
$$
The factor system defined in {ref}`space_group_irreps` satisfies this condition.

## Obtaining all irreps from (projective) regular representation

Refs. {cite}`THOMAS201776,doi:10.1137/090779966,Po2017`

For a finite group $G$ and its (projective) regular representation $\\Gamma^{(\\mathrm{reg})}$,
$$
\\tilde{\\mathbf{H}}
= \\sum\_{R \\in G} \\mathbf{\\Gamma}^{(\\mathrm{reg})}(R) \\mathbf{H} \\mathbf{\\Gamma}^{(\\mathrm{reg})}(R^{-1})
$$
is also Hermite and commute with $\\Gamma^{(\\mathrm{reg})}$ [^footnote0], where $\\mathbf{H}$ is any Hermite matrix.
Consider the spectral decomposition of $\\tilde{\\mathbf{H}}$,
$$
\\tilde{\\mathbf{H}} = \\sum\_{\\lambda} \\lambda \\sum\_{n=1}^{ d\_{\\lambda} } \\mathbf{v}_{\\lambda n} \\mathbf{v}_{\\lambda n}^{\\dagger} \\
\\lambda \\in \\mathbb{R}, \\mathbf{v}_{\\lambda n}^{\\dagger} \\mathbf{v}_{\\lambda' n'} = \\delta\_{\\lambda \\lambda'} \\delta\_{n n'}.
$$
Then, the regular representation is block-diagonalized with $\\mathbf{V} = ( \\mathbf{v}_{\\lambda_{1} 1} \\dots \\mathbf{v}_{\\lambda_{1} d\_{\\lambda\_{1}}} \\dots )$,
$$
\\left[ \\mathbf{V}^{\\dagger} \\mathbf{\\Gamma}^{(\\mathrm{reg})}(R) \\mathbf{V} \\right]\_{\\lambda n, \\lambda' n'} (\\lambda - \\lambda') = 0
$$

\[^footnote0\]: The derivation is as follows:
`{math}       \mathbf{\Delta}^{(\mathrm{reg})}(S) \mathbf{U} \mathbf{\Delta}^{(\mathrm{reg})}(S)^{-1}       &= \mathbf{\Delta}^{(\mathrm{reg})}(S) \sum_{P \in H} \mathbf{\Delta}^{(\mathrm{reg})}(P) \mathbf{\Delta}^{(\mathrm{reg})}(P)^{-1} \mathbf{\Delta}^{(\mathrm{reg})}(S)^{-1} \\       &= \sum_{P \in H} \mu(S, P) \mathbf{\Delta}^{(\mathrm{reg})}(SP) \left( \mu(S, P) \mathbf{\Delta}^{(\mathrm{reg})}(SP) \right)^{-1} \\       &= \sum_{P \in H} \mathbf{\Delta}^{(\mathrm{reg})}(P) \mathbf{\Delta}^{(\mathrm{reg})}(P)^{-1} \\       &= \mathbf{U}     `
Be careful that $\\mathbf{\\Delta}^{(\\mathrm{reg})}(S^{-1}) = \\mu(S, S^{-1}) \\mu(E, E) \\mathbf{\\Delta}^{(\\mathrm{reg})}(S)^{-1}$.

$\\mathbf{V}_{\\lambda} = ( \\mathbf{v}_{\\lambda 1} \\dots \\mathbf{v}_{\\lambda_{1} d\_{\\lambda}} ) \\in \\mathbb{C}^{ |G| \\times d\_{\\lambda} }$ forms Irrep,
$$
\\mathbf{\\Gamma}^{(\\lambda)}(R)
= \\mathbf{V}_{\\lambda}^{\\dagger} \\mathbf{\\Gamma}^{(\\mathrm{reg})}(R) \\mathbf{V}_{\\lambda}
\\quad
\\in \\mathbb{C}^{ d\_{\\lambda} \\times d\_{\\lambda} }
$$
The regular representation contains all Irreps of the group up to unitary transformation as
$$
\\Gamma^{(\\mathrm{reg})} = \\sum\_{\\alpha} d\_{\\alpha} \\Gamma^{(\\alpha)},
$$
where $d\_{\\alpha}$ is dimension of Irrep labeled as $\\alpha$, and the summation is taken over all Irreps up to unitary transformation.
Thus, this procedure produces all Irreps exhaustively.
For an infinite group, the above procedure also holds replacing regular representation into projective regular representation.

Also, we can check the obtained Irreps are enough by checking the following equality

$$
\\sum\_{\\alpha} d\_{\\alpha}^{2} &= |G| \\
\\sum\_{ g \\in G } | \\chi^{(\\alpha)}(g) |^2 &= |G| \\
\\sum\_{g \\in G} \\chi^{(\\alpha)}(g)^{\\ast} \\chi^{(\\beta)}(g) &= |G| \\delta\_{\\alpha, \\beta},
$$ (eqn:all_irreps_dim_sum)

where $\\chi^{(\\alpha)}$ is character of irrep $\\Gamma^{(\\alpha)}$.
Note that Eqs. {eq}`eqn:all_irreps_dim_sum` also holds for projective representations.

## Working example

### Crystallographic point group

#### $C\_{3v}$ associated with $P3m1$ (No. 156)

The matrix representation of this crystallographic point group is

$$
\\mathcal{P}
&=
\\left{
g_0=e,
g_1,
g_2=g_1^{-1},
g_3,
g_4=g_1^{-1} * g_3,
g_5=g_1 * g_3
\\right} \\
\\quad \\mbox{where}\\quad
g_1
&=
\\begin{pmatrix}
0 & -1 & 0 \\
1 & -1 & 0 \\
0 & 0 & 1 \\
\\end{pmatrix} \\
g_3
&=
\\begin{pmatrix}
0 & -1 & 0\\
-1 & 0 & 0\\
0 & 0 & 1\\
\\end{pmatrix}.
$$

$$
\\Gamma^{(\\mathrm{reg})} = \\Gamma^{(A\_{1})} + \\Gamma^{(A\_{2})} + 2\\Gamma^{(E)}
$$

### Space group

#### $P4\_{2}/mnm$ (No. 136) at $X=(0\\frac{1}{2}0)$

The little co-group is
$$
\\overline{\\mathcal{G}}^{\\mathbf{k}}
&= \\left{
g\_{0} = e,
g\_{1},
g\_{2},
g\_{3} = g\_{1} g\_{2},
g\_{4},
g\_{5} = g\_{1} g\_{4},
g\_{6} = g\_{2} g\_{4},
g\_{7} = g\_{1} g\_{2} g\_{4}
\\right}
\\cong mmm, \\
\\mathrm{where} \\quad
g\_{1} &= \\mathrm{diag}(-1, -1, -1) \\
g\_{2} &= \\mathrm{diag}(-1, -1, 1) \\
g\_{4} &= \\mathrm{diag}(1, -1, -1).
$$

$$
\\Gamma^{(\\mathrm{reg})} = 2 \\Gamma^{(X\_{1})} + 2 \\Gamma^{(X\_{2})}
$$

#### $\\mathcal{G} = Ia\\overline{3}d$ (No. 230) at $H=(\\frac{1}{2}\\overline{\\frac{1}{2}}\\frac{1}{2})\_{\\mathrm{primitive}}$

(corresponding to $G^{4}\_{96}$ in Ref. {cite}`Bradley2009-ze`)

$$
\\overline{\\mathcal{G}}^{\\mathbf{k}}| \\cong m\\overline{3}m,
\\quad
\\left| \\overline{\\mathcal{G}}^{\\mathbf{k}} \\right| = 48
$$

$$
\\Gamma^{(\\mathrm{reg})} = 2\\Gamma^{(H\_{1})} + 2\\Gamma^{(H\_{2})} + 6\\Gamma^{(H\_{3})}
$$

## References

```{bibliography}
:filter: docname in docnames
```
