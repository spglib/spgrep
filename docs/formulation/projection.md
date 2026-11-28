# Constructing symmetry-adapted basis

Understand how projection extracts basis functions that transform as a chosen irrep.

!!! info "Before you start"

    Read [space-group irreps](irreps/spacegroup_irreps.md) and the [array conventions](../conventions.md). You need a representation of the physical objects being projected.
    For executable applications, use the [tensor](../examples/symmetry_adapted_tensor.md) or [phonon](../examples/lattice_vibration.md) tutorial.

A typical procedure to use irreps is as follows:

1. Define an action of symmetry operations of a space group $\mathcal{G}$ on your interested objects

1. Fourier-transform your selected basis $\{ \phi^{(\mathbf{k})}_{i} \}_{\mathbf{k}, i}$ such that

    $$
    \begin{aligned}
    g \phi^{(\mathbf{k})}_{j} &= \sum_{i} \phi^{(\mathbf{k})}_{i} \Gamma^{(\mathbf{k})}_{ij}
    \quad (g \in \mathcal{G}^{\mathbf{k}}) \\
    \mathbf{\Gamma}^{(\mathbf{k})}((\mathbf{E}, \mathbf{t})) &= e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}
    \quad ( (\mathbf{E}, \mathbf{t}) \in \mathcal{G}^{\mathbf{k}}) \\
    \end{aligned}
    $$

1. Compute unitary small representations $\Gamma^{ (\mathbf{k}, \alpha) }$ of little group $\mathcal{G}^{\mathbf{k}}$ by [`spgrep.get_spacegroup_irreps_from_primitive_symmetry`][spgrep.get_spacegroup_irreps_from_primitive_symmetry] in primitive cell

1. Apply projection operator by [`spgrep.rep.irreps.project_to_irrep`][spgrep.rep.irreps.project_to_irrep]

## Projection operator

Let $\Delta^{(\alpha)}$ be **unitary** projective irrep of group $G$ with $\mu(E, E)=1$.
The projection operator can be defined as the same form with the ordinary representation [Altmann, 1989](../references.md#https://doi.org/10.1002/qua.560350309) [^proj_derivation],

$$
\begin{aligned}
P^{(\alpha)}_{ij}
    := \frac{ d_{\alpha} }{|G|} \sum_{ R \in G } \Delta^{ (\alpha) }(R)_{ij}^{\ast} R.
\end{aligned}
$$

Let $\phi^{(\alpha,j)}_{i} := P^{(\alpha)}_{ij} \phi$.
Basis vectors $\{ \phi^{(\alpha,j)}_{i} \}_{i}$ forms $\Delta^{ (\alpha) }$.

$$
\begin{aligned}
P^{(\alpha)}_{ij} P^{(\alpha)}_{i'j'}
    &= \frac{ d_{\alpha}^{2} }{|G|^{2}}
       \sum_{ R, R' \in G } \Delta^{ (\alpha) }(R)_{ij}^{\ast} \Delta^{ (\alpha) }(R)_{i'j'}^{\ast} RR' \\
    &= \frac{ d_{\alpha}^{2} }{|G|^{2}}
       \sum_{ R, T \in G } \sum_{l}
            \frac{ \mu(R, R^{1})\mu(R, R^{-1}T) }{ \mu(R^{-1}, T) }
            \Delta^{ (\alpha) }(R)_{ij}^{\ast} \Delta^{ (\alpha) }(R)_{li'} \Delta^{ (\alpha) }(T)_{lj'}^{\ast} T \\
    &= \delta_{ji'} \frac{ d_{\alpha} }{|G|}
       \sum_{ T \in G } \Delta^{ (\alpha) }(T)_{ij'}^{\ast} T
       \quad (\because \mu(E, T) = 1) \\
    &= \delta_{ji'} P^{(\alpha)}_{ij'} \\
P^{(\alpha) \dagger}_{ij}
    &= P^{(\alpha)}_{ji} \\
\end{aligned}
$$

Basis vectors $\{ \phi^{(\alpha,j)}_{i} \}_{i}$ are mutually orthogonal:

$$
\begin{aligned}
( \phi^{(\alpha,j)}_{i}, \phi^{(\alpha,j)}_{i'} )
    &= (P^{(\alpha)}_{ij}\phi, P^{(\alpha)}_{i'j}\phi) \\
    &= (\phi, P^{(\alpha)}_{ji}P^{(\alpha)}_{i'j}\phi) \\
    &= \delta_{ii'}(\phi, P^{(\alpha)}_{jj}\phi) \\
    &= \delta_{ii'}(\phi, \phi^{(\alpha, j)}_{j}).
\end{aligned}
$$

[^proj_derivation]: $\{ \phi^{ (\alpha,j,n) }_{i} \}_{i=1}^{d_{\alpha}}$ forms basis for $\Delta^{ (\alpha) }$ if $\phi^{ (\alpha,j,n) }_{i} \neq 0$.
    The derivation is as follows:

    $$
    \begin{aligned}
    S \phi^{ (\alpha,j,n) }_{i}
        &= \frac{ d_{\alpha} }{|G|} \sum_{ R \in G } \Delta^{ (\alpha) }(R)_{ij}^{\ast} R \phi \\
        &= \frac{ d_{\alpha} }{|G|} \sum_{ T \in G } \Delta^{ (\alpha) }(S^{-1}T)_{ij}^{\ast} T \phi \\
        &= \frac{ d_{\alpha} }{|G|} \sum_{ T \in G } \sum_{l} \frac{\mu(T^{-1}, S) \mu(T, T^{-1})}{\mu(S^{-1}T, T^{-1}S)} \Delta^{ (\alpha) }(S)_{li} \Delta^{ (\alpha) }(T)_{lj}^{\ast} T \phi \\
        &= \sum_{l} \phi^{ (\alpha,j,n) }_{l} \Delta^{ (\alpha) }(S)_{li}.
    \end{aligned}
    $$

    Here we assume the factor system is chosen as $\mu(E, E) = 1$.
    Then, we use $\mu(S, S^{-1}) = \mu(S^{-1}, S)$.

## Next steps

Apply the projection in the [tensor tutorial](../examples/symmetry_adapted_tensor.md) or the [phonon tutorial](../examples/lattice_vibration.md). The [representation API](../api/api_rep.md#spgrep.rep.irreps.project_to_irrep) documents `project_to_irrep`.
