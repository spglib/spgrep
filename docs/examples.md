# Examples of applications of irreps

A typical procedure to use irreps is as follows:

1. Define an action of symmetry operations of a space group {math}`\mathcal{G}` on your interested objects
2. Fourier-transform your selected basis {math}`\{ \phi^(\mathbf{k})_{i} \}_{\mathbf{k}, i}` such that
    ```{math}
        g \phi^{(\mathbf{k})}_{j} &= \sum_{i} \phi^{(\mathbf{k})}_{i} \Gamma^{(\mathbf{k})}_{ij} 
        \quad (g \in \mathcal{G}) \\
        \Gamma^{(\mathbf{k})}((E, \mathbf{t})) &= e^{ -i\mathbf{k}\cdot\mathbf{t} } \mathbf{1}
        \quad ( (E, \mathbf{t}) \in \mathcal{G}) \\
    ```
3. Compute irreps {math}`\overline{\Gamma}^{ (\mathbf{k}, \alpha) }` of little co-group {math}`\overline{\mathcal{G}}^{\mathbf{k}}` by {func}`spgrep.irreps.get_irreps_from_solvable_group_chain` or {func}`spgrep.irreps.get_irreps_from_regular`
4. Compute small representations {math}`\Gamma^{ (\mathbf{k}, \alpha) }` of little group {math}`\mathcal{G}^{\mathbf{k}}` by
    ```{math}
        \Gamma^{ (\mathbf{k}, \alpha) }( (\mathbf{R}, \mathbf{v}) )
            = e^{ -i\mathbf{k}\cdot\mathbf{v} } \overline{\Gamma}^{ (\mathbf{k}, \alpha) } ( (\mathbf{R}, \mathbf{v}) )
    ```
