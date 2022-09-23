# Constructing symmetry-adapted basis

A typical procedure to use irreps is as follows:

1. Define an action of symmetry operations of a space group {math}`\mathcal{G}` on your interested objects
1. Fourier-transform your selected basis {math}`\{ \phi^{(\mathbf{k})}_{i} \}_{\mathbf{k}, i}` such that
    ```{math}
        g \phi^{(\mathbf{k})}_{j} &= \sum_{i} \phi^{(\mathbf{k})}_{i} \Gamma^{(\mathbf{k})}_{ij} 
        \quad (g \in \mathcal{G}^{\mathbf{k}}) \\
        \Gamma^{(\mathbf{k})}((E, \mathbf{t})) &= e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}
        \quad ( (E, \mathbf{t}) \in \mathcal{G}^{\mathbf{k}}) \\
    ```
1. Compute unitary small representations {math}`\Gamma^{ (\mathbf{k}, \alpha) }` of little group {math}`\mathcal{G}^{\mathbf{k}}` by {func}`spgrep.get_spacegroup_irreps_from_primitive_symmetry` in primitive cell
1. Apply projection operator by {func}`spgrep.representation.project_to_irrep`

## Projection operator

Let {math}`\Delta^{(\alpha)}` be **unitary** projective irrep of group {math}`G` with {math}`\mu(E, E)=1`.
The projection operator can be defined as the same form with the ordinary representation {cite}`https://doi.org/10.1002/qua.560350309` [^proj_derivation],
```{math}
    P^{(\alpha)}_{ij}
        := \frac{ d_{\alpha} }{|G|} \sum_{ R \in G } \Delta^{ (\alpha) }(R)_{ij}^{\ast} R.
```
Let {math}`\phi^{(\alpha,j)}_{i} := P^{(\alpha)}_{ij} \phi`.
Basis vectors {math}`\{ \phi^{(\alpha,j)}_{i} \}_{i}` forms {math}`\Delta^{ (\alpha) }`.

```{math}
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
```

Basis vectors {math}`\{ \phi^{(\alpha,j)}_{i} \}_{i}` are mutually orthogonal:
```{math}
    ( \phi^{(\alpha,j)}_{i}, \phi^{(\alpha,j)}_{i'} )
        &= (P^{(\alpha)}_{ij}\phi, P^{(\alpha)}_{i'j}\phi) \\
        &= (\phi, P^{(\alpha)}_{ji}P^{(\alpha)}_{i'j}\phi) \\
        &= \delta_{ii'}(\phi, P^{(\alpha)}_{jj}\phi) \\
        &= \delta_{ii'}(\phi, \phi^{(\alpha, j)}_{j}) \\
```

[^proj_derivation]: {math}`\{ \phi^{ (\alpha,j,n) }_{i} \}_{i=1}^{d_{\alpha}}` forms basis for {math}`\Delta^{ (\alpha) }` if {math}`\phi^{ (\alpha,j,n) }_{i} \neq 0`.
    The derivation is as follows:
    ```{math}
        S \phi^{ (\alpha,j,n) }_{i}
            &= \frac{ d_{\alpha} }{|G|} \sum_{ R \in G } \Delta^{ (\alpha) }(R)_{ij}^{\ast} R \phi \\
            &= \frac{ d_{\alpha} }{|G|} \sum_{ T \in G } \Delta^{ (\alpha) }(S^{-1}T)_{ij}^{\ast} T \phi \\
            &= \frac{ d_{\alpha} }{|G|} \sum_{ T \in G } \sum_{l} \frac{\mu(T^{-1}, S) \mu(T, T^{-1})}{\mu(S^{-1}T, T^{-1}S)} \Delta^{ (\alpha) }(S)_{li} \Delta^{ (\alpha) }(T)_{lj}^{\ast} T \phi \\
            &= \sum_{l} \phi^{ (\alpha,j,n) }_{l} \Delta^{ (\alpha) }(S)_{li}.
    ```
    Here we assume the factor system is chosen as {math}`\mu(E, E) = 1`.
    Then, we use {math}`\mu(S, S^{-1}) = \mu(S^{-1}, S)`.

## Symmetry-adapted tensor with intrinsic symmetry

{cite}`el2008symmetry`

Consider vector space {math}`V` and its symmetry adapted basis {math}`\{ \mathbf{f}^{\alpha m}_{i} \}` w.r.t. group {math}`G`.
```{math}
    V &= \bigoplus_{\alpha} \bigoplus_{m=1}^{m_{\alpha}} V^{\alpha m} \\
    V^{\alpha m} &= \bigoplus_{i=1}^{d_{\alpha}} K \mathbf{f}^{\alpha m}_{i} \\
    g \mathbf{f}^{\alpha m}_{j}
        &= \sum_{i=1} ^{d_{\alpha}} \mathbf{f}^{\alpha m}_{i} \Gamma^{\alpha}_{ij}(g)
        \quad (g \in G, j = 1, \dots, d_{\alpha}),
```
where {math}`K` is {math}`\mathbb{C}` or {math}`\mathbb{R}`, and {math}`\Gamma^{\alpha}` is irrep over {math}`K`.

Action of {math}`G` on rank-{math}`p` tensor {math}`\mathsf{T}: V^{\ast \otimes p}` is defined as
```{math}
    (g \mathsf{T})(\mathbf{v}_{1}, \dots, \mathbf{v}_{p})
        := \mathsf{T}(g^{-1} \mathbf{v}_{1}, \dots, g^{-1} \mathbf{v}_{p})
        \quad (g \in G).
```
We also consider intrinsic symmetry {math}`\Sigma` of {math}`\mathsf{T}` as [^check_action]
```{math}
    (\sigma \mathsf{T})(\mathbf{v}_{1}, \dots, \mathbf{v}_{p})
        := \mathsf{T}(\mathbf{v}_{\sigma^{-1}(1)}, \dots, \mathbf{v}_{\sigma^{-1}(p)})
        \quad (\sigma \in \Sigma).
```

[^check_action]: These definitions follow the condition of left actions:
    ```{math}
        (\sigma (\sigma' \mathsf{T}))(\mathbf{v}_{1}, \dots, \mathbf{v}_{p})
            &= (\sigma' \mathsf{T})(\mathbf{v}_{\sigma^{-1}(1)}, \dots, \mathbf{v}_{\sigma^{-1}(p)}) \\
            &= \mathsf{T}(\mathbf{v}_{\sigma'^{-1}(\sigma^{-1}(1))}, \dots, \mathbf{v}_{\sigma'^{-1}(\sigma^{-1}(p))}) \\
            &= \mathsf{T}(\mathbf{v}_{\sigma\sigma'(1)}, \dots, \mathbf{v}_{\sigma\sigma'(p)}) \\
            &= ((\sigma \sigma') \mathsf{T})(\mathbf{v}_{1}, \dots, \mathbf{v}_{p}) \\
        \therefore \sigma (\sigma' \mathsf{T}) &= (\sigma \sigma') \mathsf{T}
    ```

## References

```{bibliography}
:filter: docname in docnames
```
