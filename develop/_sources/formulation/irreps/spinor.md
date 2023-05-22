(spin_representation)=
# Spin representation

A spin representation is a projective representation with a factor system for spinor wave functions.
Here, we give the convention to choose the factor system for spinor in spgrep.

## Factor system for spinor

### Symmetry operation of the first kind

The map $\mathrm{SO}(3) \ni \mathbf{R}_{ \theta\hat{\mathbf{n}} } \mapsto \mathbf{U} ( \mathbf{R}_{ \theta\hat{\mathbf{n}} } ) := \exp \left( -\frac{i}{2}\theta \hat{\mathbf{n}} \cdot \mathbf{\sigma} \right) \in \mathrm{SU}(2)$ is not surjective, where $\mathbf{\sigma}_{i} \, (i=x,y,z)$ are Pauli matrices
$$
\mathbf{\sigma}_{x}
    = \begin{pmatrix}
        0 & 1 \\
        1 & 0 \\
    \end{pmatrix},
\mathbf{\sigma}_{y}
    = \begin{pmatrix}
        0 & -i \\
        i & 0 \\
    \end{pmatrix},
\mathbf{\sigma}_{z}
    = \begin{pmatrix}
        1 & 0 \\
        0 & -1 \\
    \end{pmatrix}.
$$
In fact, $\mathbf{R}_{\theta\hat{\mathbf{n}}}$ and $\mathbf{R}_{2\pi - \theta, -\hat{\mathbf{n}}}$ represents the identical rotation.
However they gives different unitary matrices with opposite signs.
We choose either of the unitary matices for each rotations as convention.
In particular, we choose $\theta=0$ for identity $\mathbf{E}$.
We define a factor system from the ambiguity as
$$
\mathbf{U}(\mathbf{R}) \mathbf{U}(\mathbf{R}')
    =: z(\mathbf{R}, \mathbf{R}') \mathbf{U}(\mathbf{R}\mathbf{R}')
    \quad (\mathbf{R}, \mathbf{R}' \in \mathrm{SO}(3)).
$$
In our convention, $z(\mathbf{E}, \mathbf{R}) = z(\mathbf{R}, \mathbf{E}) = 1 \,(\forall \mathbf{R} \in \mathrm{SO}(3))$.
Also, this representation matrix adapts Condon-Shortley phase.

We define the action of a symmetry operation of the first kind on spinor as
$$
\left[ (\mathbf{R}, \mathbf{v}) \mathbf{\Psi} \right] (\mathbf{r})
    &:= \mathbf{U}( \mathbf{R} ) \mathbf{\Psi}( (\mathbf{R}, \mathbf{v})^{-1} \mathbf{r}) \\
\Rightarrow \left[ g_{1} \left[ g_{2} \mathbf{\Psi} \right] \right](\mathbf{r})
    &= \mathbf{U}(\mathbf{p}_{g_{1}}) \left[ g_{2} \mathbf{\Psi} \right](g_{1}^{-1}\mathbf{r}) \\
    &= z(\mathbf{p}_{g_{1}}, \mathbf{p}_{g_{2}}) \left[ (g_{1}g_{2}) \mathbf{\Psi} \right](\mathbf{r}).
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
        \quad (\mathbf{t} \in L_{\mathcal{T}}).
$$

A transformed Bloch function $(\mathbf{R}, \mathbf{v}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})$ is a Bloch function with $\mathbf{Rk}$:
$$
(\mathbf{E}, \mathbf{t}) (\mathbf{R}, \mathbf{v}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})
    &= (\mathbf{R}, \mathbf{v}) (\mathbf{E}, \mathbf{R}^{-1}\mathbf{t}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})
        \quad (\because z(\mathbf{E}, \mathbf{R}) = z(\mathbf{R}, \mathbf{E}) = 1 ) \\
    &= \exp \left( -i \mathbf{k} \cdot \mathbf{R}^{-1} \mathbf{t} \right) (\mathbf{R}, \mathbf{v}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r}) \\
    &= \exp \left( -i \mathbf{Rk} \cdot \mathbf{t} \right) (\mathbf{R}, \mathbf{v}) \mathbf{\Psi}_{\mathbf{k}}(\mathbf{r})
        \quad (\because \mathbf{R} \in \mathrm{SO}(3) ).
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
  \mathbf{\Gamma}^{\mathbf{k}\alpha}((\mathbf{R}, \mathbf{v}))
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
    \quad (\mathbf{R}, \mathbf{R}' \in \mathrm{SO}(3)).
$$

## Convention of rotations for spinor

A rotation $\mathbf{R}_{ \theta\hat{\mathbf{n}} } \in \mathrm{SO}(3)$ can be written with angular momentum operators:
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

The corresponding unitary rotation is explicitly written as
$$
\mathbf{U} ( \mathbf{R}_{ \theta\hat{\mathbf{n}} } )
    &= \exp \left( -\frac{i}{2}\theta \hat{\mathbf{n}} \cdot \mathbf{\sigma} \right) \\
    &= \mathbf{1} \cos \frac{\theta}{2} - i (\hat{\mathbf{n}} \cdot \mathbf{\sigma}) \sin \frac{\theta}{2} \\
    &=
        \begin{pmatrix}
            \cos \frac{\theta}{2} - i n_{z} \sin \frac{\theta}{2} & (-i n_{x} - n_{y}) \sin \frac{\theta}{2} \\
            (-i n_{x} + n_{y}) \sin \frac{\theta}{2} & \cos \frac{\theta}{2} + i n_{z} \sin \frac{\theta}{2} \\
        \end{pmatrix}.
$$

## References

```{bibliography}
:filter: docname in docnames
```
