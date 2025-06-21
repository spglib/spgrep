# Memo

- Spin-orbit coupling in VASP: {cite}`PhysRevB.62.11556`, {cite}`PhysRevB.93.224425`
- $SU(2)$ and $SO(3)$: {cite}`altmann2005rotations`, {cite}`so3su2`

## Bloch sphere

Consider a spinor $\\psi\_{\\uparrow}(\\mathbf{r})\\ket{\\uparrow} + \\psi\_{\\downarrow}(\\mathbf{r})\\ket{\\downarrow}$ with $|\\psi\_{\\uparrow}(\\mathbf{r})|^{2} + |\\psi\_{\\downarrow}(\\mathbf{r})|^{2} = 1$.

The following construction maps the spinor to a point $\\left( m\_{x}(\\mathbf{r}), m\_{y}(\\mathbf{r}), m\_{z}(\\mathbf{r}) \\right) \\in S^{2}$:
$$
m\_{x}(\\mathbf{r})
&:=
\\begin{pmatrix} \\psi\_{\\uparrow}(\\mathbf{r})^{\\ast} & \\psi\_{\\downarrow}(\\mathbf{r})^{\\ast} \\end{pmatrix}
\\mathbf{\\sigma}_{x}
\\begin{pmatrix} \\psi_{\\uparrow}(\\mathbf{r}) \\ \\psi\_{\\downarrow}(\\mathbf{r}) \\end{pmatrix} \\
&=
2 ,\\mathrm{Re} \\left( \\psi\_{\\uparrow}(\\mathbf{r})^{\\ast} \\psi\_{\\downarrow}(\\mathbf{r}) \\right) \\in \\mathbb{R} \\
m\_{y}(\\mathbf{r})
&:=
\\begin{pmatrix} \\psi\_{\\uparrow}(\\mathbf{r})^{\\ast} & \\psi\_{\\downarrow}(\\mathbf{r})^{\\ast} \\end{pmatrix}
\\mathbf{\\sigma}_{y}
\\begin{pmatrix} \\psi_{\\uparrow}(\\mathbf{r}) \\ \\psi\_{\\downarrow}(\\mathbf{r}) \\end{pmatrix} \\
&=
2 ,\\mathrm{Im} \\left( \\psi\_{\\uparrow}(\\mathbf{r})^{\\ast} \\psi\_{\\downarrow}(\\mathbf{r}) \\right) \\in \\mathbb{R} \\
m\_{z}(\\mathbf{r})
&:=
\\begin{pmatrix} \\psi\_{\\uparrow}(\\mathbf{r})^{\\ast} & \\psi\_{\\downarrow}(\\mathbf{r})^{\\ast} \\end{pmatrix}
\\mathbf{\\sigma}_{z}
\\begin{pmatrix} \\psi_{\\uparrow}(\\mathbf{r}) \\ \\psi\_{\\downarrow}(\\mathbf{r}) \\end{pmatrix} \\
&=
|\\psi\_{\\uparrow}(\\mathbf{r})|^{2} - |\\psi\_{\\downarrow}(\\mathbf{r})|^{2} \\in \\mathbb{R} \\
$$

$$
m\_{x}(\\mathbf{r})^{2} + m\_{y}(\\mathbf{r})^{2} + m\_{z}(\\mathbf{r})^{2} = 1
$$

## Symmetry-adapted tensor with intrinsic symmetry

Ref. {cite}`el2008symmetry`

Consider vector space $V$ and its symmetry adapted basis ${ \\mathbf{f}^{\\alpha m}_{i} }$ w.r.t. group $G$.
$$
V &= \\bigoplus_{\\alpha} \\bigoplus\_{m=1}^{m\_{\\alpha}} V^{\\alpha m} \\
V^{\\alpha m} &= \\bigoplus\_{i=1}^{d\_{\\alpha}} K \\mathbf{f}^{\\alpha m}_{i} \\
g \\mathbf{f}^{\\alpha m}_{j}
&= \\sum\_{i=1} ^{d\_{\\alpha}} \\mathbf{f}^{\\alpha m}_{i} \\Gamma^{\\alpha}_{ij}(g)
\\quad (g \\in G, j = 1, \\dots, d\_{\\alpha}),
$$
where $K$ is $\\mathbb{C}$ or $\\mathbb{R}$, and $\\Gamma^{\\alpha}$ is irrep over $K$.

Action of $G$ on rank-$p$ tensor $\\mathsf{T}: V^{\\ast \\otimes p}$ is defined as
$$
(g \\mathsf{T})(\\mathbf{v}_{1}, \\dots, \\mathbf{v}_{p})
:= \\mathsf{T}(g^{-1} \\mathbf{v}_{1}, \\dots, g^{-1} \\mathbf{v}_{p})
\\quad (g \\in G).
$$
We also consider intrinsic symmetry $\\Sigma$ of $\\mathsf{T}$ as [^check_action]
$$
(\\sigma \\mathsf{T})(\\mathbf{v}_{1}, \\dots, \\mathbf{v}_{p})
:= \\mathsf{T}(\\mathbf{v}_{\\sigma^{-1}(1)}, \\dots, \\mathbf{v}_{\\sigma^{-1}(p)})
\\quad (\\sigma \\in \\Sigma).
$$

\[^check_action\]: These definitions follow the condition of left actions:
`{math}         (\sigma (\sigma' \mathsf{T}))(\mathbf{v}_{1}, \dots, \mathbf{v}_{p})             &= (\sigma' \mathsf{T})(\mathbf{v}_{\sigma^{-1}(1)}, \dots, \mathbf{v}_{\sigma^{-1}(p)}) \\             &= \mathsf{T}(\mathbf{v}_{\sigma'^{-1}(\sigma^{-1}(1))}, \dots, \mathbf{v}_{\sigma'^{-1}(\sigma^{-1}(p))}) \\             &= \mathsf{T}(\mathbf{v}_{\sigma\sigma'(1)}, \dots, \mathbf{v}_{\sigma\sigma'(p)}) \\             &= ((\sigma \sigma') \mathsf{T})(\mathbf{v}_{1}, \dots, \mathbf{v}_{p}) \\         \therefore \sigma (\sigma' \mathsf{T}) &= (\sigma \sigma') \mathsf{T}     `

## References

```{bibliography}
:filter: docname in docnames
```
