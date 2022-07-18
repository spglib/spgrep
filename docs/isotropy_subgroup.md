# Isotropy subgroup

Consider space group {math}`\mathcal{G}` and its representation
```{math}
    g \phi_{j} = \sum_{i} \phi_{i} \Gamma(g)_{ij} \quad (g \in \mathcal{G}).
```
A subgroup {math}`\mathcal{H} (\leq \mathcal{G})` is called **isotropy subgroup** belonging {math}`\Gamma` if the subduced representation {math}`\Gamma \downarrow \mathcal{H}` is identity [^isotropy_subgroup] {cite}`Howard:ta0003,PhysRevB.43.11010,doi:10.1142/0751,Stokes:vk5013`.

When we consider a linear combination of basis functions as {math}`\mathbf{\eta} = \sum_{i} \eta_{i} \phi_{i}`, we call {math}`\mathbf{\eta} = \{ \eta_{i} \}_{i}` as **order parameters**.
Two order parameters {math}`\mathbf{\eta}` and {math}`\mathbf{\eta}'` are equivalent if some operation {math}`g \in \mathcal{G}` exists such that
```{math}
    \mathbf{\eta}' = \mathbf{\Gamma}(g) \mathbf{\eta}.
```

[^isotropy_subgroup]: Consider group {math}`G` acting on space {math}`X`.
    In general, isotropy subgroup is defined as stabilizer of point {math}`x` in {math}`X` as
    ```{math}
        G_{x} = \left\{ g \in G \mid g x = x \right\}.
    ```

## Asymmetric units

Ref. {cite}`Grosse-Kunstleve:pz5088`

## References

```{bibliography}
:filter: docname in docnames
```
