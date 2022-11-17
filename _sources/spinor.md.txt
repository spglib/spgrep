# Spin representation

- Spin-orbit coupling in VASP: {cite}`PhysRevB.62.11556`, {cite}`PhysRevB.93.224425`
- $SU(2)$ and $SO(3)$: {cite}`altmann2005rotations`, {cite}`so3su2`

## Bloch sphere

Consider a spinor $\psi_{\uparrow}(\mathbf{r})\ket{\uparrow} + \psi_{\downarrow}(\mathbf{r})\ket{\downarrow}$ with $|\psi_{\uparrow}(\mathbf{r})|^{2} + |\psi_{\downarrow}(\mathbf{r})|^{2} = 1$.

The following construction maps the spinor to a point $\left( m_{x}(\mathbf{r}), m_{y}(\mathbf{r}), m_{z}(\mathbf{r}) \right) \in S^{2}$:
$$
m_{x}(\mathbf{r})
    &:=
    \begin{pmatrix} \psi_{\uparrow}(\mathbf{r})^{\ast} & \psi_{\downarrow}(\mathbf{r})^{\ast} \end{pmatrix}
    \mathbf{\sigma}_{x}
    \begin{pmatrix} \psi_{\uparrow}(\mathbf{r}) \\ \psi_{\downarrow}(\mathbf{r}) \end{pmatrix} \\
    &=
    2 \,\mathrm{Re} \left( \psi_{\uparrow}(\mathbf{r})^{\ast} \psi_{\downarrow}(\mathbf{r}) \right) \in \mathbb{R} \\
m_{y}(\mathbf{r})
    &:=
    \begin{pmatrix} \psi_{\uparrow}(\mathbf{r})^{\ast} & \psi_{\downarrow}(\mathbf{r})^{\ast} \end{pmatrix}
    \mathbf{\sigma}_{y}
    \begin{pmatrix} \psi_{\uparrow}(\mathbf{r}) \\ \psi_{\downarrow}(\mathbf{r}) \end{pmatrix} \\
    &=
    2 \,\mathrm{Im} \left( \psi_{\uparrow}(\mathbf{r})^{\ast} \psi_{\downarrow}(\mathbf{r}) \right) \in \mathbb{R} \\
m_{z}(\mathbf{r})
    &:=
    \begin{pmatrix} \psi_{\uparrow}(\mathbf{r})^{\ast} & \psi_{\downarrow}(\mathbf{r})^{\ast} \end{pmatrix}
    \mathbf{\sigma}_{z}
    \begin{pmatrix} \psi_{\uparrow}(\mathbf{r}) \\ \psi_{\downarrow}(\mathbf{r}) \end{pmatrix} \\
    &=
    |\psi_{\uparrow}(\mathbf{r})|^{2} - |\psi_{\downarrow}(\mathbf{r})|^{2} \in \mathbb{R} \\
$$
$$
m_{x}(\mathbf{r})^{2} + m_{y}(\mathbf{r})^{2} + m_{z}(\mathbf{r})^{2} = 1.
$$

(spinor_factor_system)=
## Factor system for spinor

### Symmetry operation of the first kind

The map $SO(3) \ni \mathbf{R}_{ \theta\hat{\mathbf{n}} } \mapsto \mathbf{U} ( \mathbf{R}_{ \theta\hat{\mathbf{n}} } ) := \exp \left( -\frac{i}{2}\theta \hat{\mathbf{n}} \cdot \mathbf{\sigma} \right) \in SU(2)$ is not surjective.
In fact, $\mathbf{R}_{\theta\hat{\mathbf{n}}}$ and $\mathbf{R}_{2\pi - \theta, -\hat{\mathbf{n}}}$ represents the identical rotation.
However they gives different unitary matrices with opposite signs.
We choose either of the unitary matices for each rotations as convention.
In particular, we choose $\theta=0$ for identity $\mathbf{E}$.
We define a factor system from the ambiguity as
$$
\mathbf{U}(g) \mathbf{U}(g') = z(g, g') \mathbf{U}(gg') \quad (g, g' \in SO(3)).
$$
In our convention, $z(\mathbf{E}, g) = z(g, \mathbf{E}) = 1 \,(\forall g \in SO(3))$.
Also, this representation matrix adapts Condon-Shortley phase.

We define the action of a symmetry operation of the first kind on spinor as
$$
(\mathbf{R}_{\theta\hat{\mathbf{n}}}, \mathbf{v}) \mathbf{\Psi}(\mathbf{r})
    := \mathbf{U}( \mathbf{R}_{ \theta \hat{\mathbf{n}} } ) \mathbf{\Psi}( (\mathbf{R}_{\theta\hat{\mathbf{n}}}, \mathbf{v})^{-1} \mathbf{r}).
$$

Consider Bloch function with $\mathbf{k}$
$$
\mathbf{\Psi}_{\mathbf{k}}(\mathbf{r}) = e^{ i \mathbf{k} \cdot \mathbf{r} } \mathbf{u}_{\mathbf{k}}(\mathbf{r}),
$$
where $\mathbf{u}_{\mathbf{k}}(\mathbf{r})$ is periodic w.r.t. lattice translations $L_{\mathcal{T}}$.
The Bloch function forms an irreps of translation subgroup $\mathcal{T}$ of space group $\mathcal{G}$:
$$
(\mathbf{E}, \mathbf{t}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})
    &= \mathbf{\Psi}( (\mathbf{E}, \mathbf{t})^{-1} \mathbf{r}) \\
    &= e^{ -i \mathbf{k} \cdot \mathbf{t} } \mathbf{\Psi}(\mathbf{r})
        \quad (\mathbf{t} \in L_{\mathcal{T}})
$$

$(\mathbf{R}, \mathbf{v}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})$ is a Bloch function with $\mathbf{Rk}$:
$$
(\mathbf{E}, \mathbf{t}) (\mathbf{R}, \mathbf{v}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})
    &= (\mathbf{R}, \mathbf{v}) (\mathbf{E}, \mathbf{R}^{-1}\mathbf{t}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})
        \quad (\because z(\mathbf{E}, \mathbf{R}) = z(\mathbf{R}, \mathbf{E}) = 1 ) \\
    &= \exp \left( -i \mathbf{k} \cdot \mathbf{R}^{-1} \mathbf{t} \right) (\mathbf{R}, \mathbf{v}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r}) \\
    &= \exp \left( -i \mathbf{Rk} \cdot \mathbf{t} \right) (\mathbf{R}, \mathbf{v}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})
        \quad (\because \mathbf{R} \in SO(3) ).
$$

Let {math}`\Gamma^{\mathbf{k}\alpha}` be a projective irrep of the little group {math}`\mathcal{G}^{\mathbf{k}} = \coprod_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } (\mathbf{S}_{i}, \mathbf{w}_{i}) \mathcal{T}` and $\{ \mathbf{\Psi}^{\mathbf{k}\alpha}_{\mu} \}_{\mu=1}^{d_{\mathbf{k}\alpha}}$ form the projective irrep.
For $\mathbf{S}_{i} \mathbf{S}_{j} = \mathbf{S}_{k}$,
$$
    &\sum_{\mu} \mathbf{\Psi}^{\mathbf{k}\alpha}_{\mu} \left[ \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{i}, \mathbf{w}_{i})) \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{j}, \mathbf{w}_{j})) \right]_{\mu\nu} \\
    &=
        (\mathbf{S}_{i}, \mathbf{w}_{i}) \left( (\mathbf{S}_{j}, \mathbf{w}_{j}) \mathbf{\Psi}^{\mathbf{k}\alpha}_{\nu} \right) \\
    &=
        \frac{ z(\mathbf{S}_{i}, \mathbf{S}_{j}) }{ z(\mathbf{E}, \mathbf{S}_{k}) }
        (\mathbf{S}_{k}, \mathbf{w}_{k})
        \left( (\mathbf{E}, \mathbf{S}_{k}^{-1} (\mathbf{S}_{i}\mathbf{w}_{j} + \mathbf{w}_{i} - \mathbf{w}_{k} ) ) \mathbf{\Psi}^{\mathbf{k}\alpha}_{\nu} \right) \\
    &=
        \sum_{\mu} \mathbf{\Psi}^{\mathbf{k}\alpha}_{\mu}
        z(\mathbf{S}_{i}, \mathbf{S}_{j})
        e^{-i\mathbf{k}\cdot \left( \mathbf{S}_{i}\mathbf{w}_{j} + \mathbf{w}_{i} - \mathbf{w}_{k} \right) } 
        \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{S}_{k}, \mathbf{w}_{k}))
        \quad (\because z(\mathbf{E}, \mathbf{S}_{k}) = 1, \mathbf{S}_{k} \in \overline{\mathcal{G}}^{\mathbf{k}}).
$$

Then, Simplifying by
$$
  \Gamma^{\mathbf{k}\alpha}((\mathbf{R}, \mathbf{v}))
    &=: e^{ -i \mathbf{k} \cdot \mathbf{v} } \mathbf{D}^{\mathbf{k}\alpha}((\mathbf{R}, \mathbf{v})) \\
  \mathbf{D}^{\mathbf{k}\alpha}((\mathbf{S}_{i}, \mathbf{w}_{i})) \mathbf{D}^{\mathbf{k}\alpha}((\mathbf{S}_{j}, \mathbf{w}_{j}))
    &= z(\mathbf{S}_{i}, \mathbf{S}_{j}) e^{ -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} } \mathbf{D}^{\mathbf{k}\alpha}((\mathbf{S}_{k}, \mathbf{w}_{k})),
$$
where $\mathbf{g}_{i} = \mathbf{S}_{i}^{-1} \mathbf{k} - \mathbf{k}$.
We can enumerate projective irreps for spinor from the factor system $\mu_{\mathbf{k}}(S_{i}, S_{j}) := z(\mathbf{S}_{i}, \mathbf{S}_{j}) \exp \left( -i \mathbf{g}_{i} \cdot \mathbf{w}_{j} \right)$.

### Symmetry operation of the second kind

We **assume** the inversion $\mathbf{I}$ acts on spin functions as
$$
\mathbf{I} \ket{\uparrow} &= \ket{\uparrow} \\
\mathbf{I} \ket{\downarrow} &= \ket{\downarrow},
$$
which is known as Pauli gauge {cite}`altmann2005rotations`.
In this convention, we can choose factor system as follows:
$$
z(\mathbf{I}, \mathbf{I}) &= 1 \\
z(\mathbf{R}, \mathbf{I}\mathbf{R}') &= z(\mathbf{IR}, \mathbf{R}') = z(\mathbf{IR}, \mathbf{I}\mathbf{R}') = z(\mathbf{R}, \mathbf{R}')
    \quad (\mathbf{R}, \mathbf{R}' \in SO(3)).
$$

### Convention of rotations for spinor in spgrep

A rotation $\mathbf{R}_{ \theta\hat{\mathbf{n}} } \in SO(3)$ can be written with angular momentum operators:
$$
\mathbf{R}_{ \theta\hat{\mathbf{n}} }
    &= \exp \left( -i \theta \hat{\mathbf{n}} \cdot \mathbf{L} \right) \\
    &= \begin{pmatrix}
        \cos \theta + n_{x}^{2} (1 - \cos \theta)         & n_{x} n_{y} (1 - \cos \theta) - n_{z} \sin \theta & n_{x} n_{z} (1 - \cos \theta) + n_{y} \sin \theta \\
        n_{y} n_{x} (1 - \cos \theta) + n_{z} \sin \theta & \cos \theta + n_{y}^{2} (1 - \cos \theta)         & n_{y} n_{z} (1 - \cos \theta) - n_{x} \sin \theta \\
        n_{z} n_{x} (1 - \cos \theta) - n_{y} \sin \theta & n_{z} n_{y} (1 - \cos \theta) + n_{x} \sin \theta & \cos \theta + n_{z}^{2} (1 - \cos \theta) \\
    \end{pmatrix} \\
\mathbf{L}_{x}
    &:= -i \begin{pmatrix}
        0 & 0 & 0 \\
        0 & 0 & 1 \\
        0 & -1 & 0 \\
    \end{pmatrix} \\
\mathbf{L}_{y}
    &:= -i \begin{pmatrix}
        0 & 0 & -1 \\
        0 & 0 & 0 \\
        1 & 0 & 0 \\
    \end{pmatrix} \\
\mathbf{L}_{z}
    &:= -i \begin{pmatrix}
        0 & 1 & 0 \\
        -1 & 0 & 0 \\
        0 & 0 & 0 \\
    \end{pmatrix} \\
$$

The axis and angle of $\mathbf{R}_{ \theta\hat{\mathbf{n}} }$ are interpret by the following relations:
$$
\mathbf{R}_{f}
    &:= \mathbf{A}^{-1} \mathbf{R}_{ \theta\hat{\mathbf{n}} } \mathbf{A}
        \quad \in \mathrm{SL}(3) \\
\mathrm{tr} \mathbf{R}_{ \theta\hat{\mathbf{n}} }
    &= 2 \cos \theta + 1 \\
$$

$$
\begin{pmatrix}
    [R_{ \theta\hat{\mathbf{n}} }]_{23} - [R_{ \theta\hat{\mathbf{n}} }]_{32} \\
    [R_{ \theta\hat{\mathbf{n}} }]_{31} - [R_{ \theta\hat{\mathbf{n}} }]_{13} \\
    [R_{ \theta\hat{\mathbf{n}} }]_{12} - [R_{ \theta\hat{\mathbf{n}} }]_{21} \\
\end{pmatrix}
= -2 (\sin \theta) \hat{\mathbf{n}}
$$


Unitary rotation

$$
\mathbf{U} ( \mathbf{R}_{ \theta\hat{\mathbf{n}} } )
    &= \exp \left( -\frac{i}{2}\theta \hat{\mathbf{n}} \cdot \mathbf{\sigma} \right) \\
    &= \mathbf{1} \cos \frac{\theta}{2} - i (\hat{\mathbf{n}} \cdot \mathbf{\sigma}) \sin \frac{\theta}{2} \\
    &=
        \begin{pmatrix}
            \cos \frac{\theta}{2} - i n_{z} \sin \frac{\theta}{2} & (-i n_{x} - n_{y}) \sin \frac{\theta}{2} \\
            (-i n_{x} + n_{y}) \sin \frac{\theta}{2} & \cos \frac{\theta}{2} + i n_{z} \sin \frac{\theta}{2} \\
        \end{pmatrix}
$$

## References

```{bibliography}
:filter: docname in docnames
```
