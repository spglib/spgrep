# Lattice vibration

Refs. {cite}`birman1974quick,RevModPhys.40.1`

## Harmonic phonon

Hamiltonian and commutation relations:
```{math}
  H = \sum_{ l \kappa \mu } \frac{ p_{\mu}(l\kappa)^{2} }{ 2 M_{\kappa} }
    + \frac{1}{2} \sum_{ l \kappa \mu } \sum_{ l' \kappa' \mu' } \Phi_{ \mu \mu' }(l\kappa, l'\kappa') u_{\mu}(l\kappa) u_{\mu'}(l'\kappa')
```

Dynamical matrix [^dynamical_matrix]
```{math}
  D_{\mu \mu'}(\kappa\kappa'; \mathbf{q})
    = \frac{1}{ \sqrt{ M_{\kappa} M_{\kappa'} } } \sum_{l'} \Phi_{ \mu \mu' }(0\kappa, l'\kappa') e^{i \mathbf{q} \cdot ( \mathbf{r}(l'\kappa') - \mathbf{r}(0\kappa) )}
```
We denote {math}`[\mathbf{D}(\mathbf{q})]_{\kappa\mu, \kappa'\mu'} = D_{\mu \mu'}(\kappa\kappa'; \mathbf{q})`, then
```{math}
  \mathbf{D}(-\mathbf{q}) &= \mathbf{D}(\mathbf{q})^{\ast} \\
  \mathbf{D}(\mathbf{q})^{\dagger} &= \mathbf{D}(\mathbf{q}) \quad \mbox{(Hermite)}.
```
Here we use {math}`\Phi_{ \mu' \mu }(0\kappa', -l\kappa) = \Phi_{ \mu \mu' }(0\kappa, l'\kappa')`.

Let normalized eigenvector of {math}`\mathbf{D}(\mathbf{q})` as {math}`[\mathbf{e}(\mathbf{q}\nu)]_{\kappa\mu} = e_{\mu}(\kappa; \mathbf{q}\nu)` with
```{math}
  \mathbf{D}(\mathbf{q}) \mathbf{e}(\mathbf{q}\nu)
    &= \omega_{\mathbf{q}\nu}^{2} \mathbf{e}(\mathbf{q}\nu)
    \quad (\nu = 1, \dots, 3N) \\
  \sum_{ \kappa\mu} \sum_{ \kappa'\mu' } e_{\mu}(\kappa; \mathbf{q}\nu)^{\ast} D_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) e_{\mu'}(\kappa'; \mathbf{q}\nu')
    &= \omega_{\mathbf{q}\nu}^{2} \delta_{\nu \nu'} \\
  \sum_{ \kappa\mu} e_{\mu}(\kappa; \mathbf{q}\nu)^{\ast} e_{\mu}(\kappa; \mathbf{q}\nu')
    &= \delta_{\nu \nu'}
    \quad \mbox{(Orthogonality)} \\
  \sum_{ \nu } e_{\mu}(\kappa; \mathbf{q}\nu)^{\ast} e_{\mu'}(\kappa'; \mathbf{q}\nu)
    &= \delta_{\kappa\kappa'} \delta_{\mu\mu'}
    \quad \mbox{(Completeness)} \\
  \omega_{-\mathbf{q}\nu}^{2} &= \omega_{\mathbf{q}\nu}^{2}.
```
We can choose as
```{math}
  e_{\mu}(\kappa; -\mathbf{q}\nu) = e_{\mu}(\kappa; \mathbf{q}\nu)^{\ast}.
```
Later we denote {math}`q = (\mathbf{q}, \nu)` and {math}`-q = (-\mathbf{q}, \nu)`.

[^dynamical_matrix]: This is the same phase convention with [phonopy](https://phonopy.github.io/phonopy/formulation.html#dynamical-matrix).
    There is the other formulation for defining dynamical matrix as
    ```{math}
    \tilde{D}_{\mu \mu'}(\kappa\kappa'; \mathbf{q})
      = \frac{1}{ \sqrt{ M_{\kappa} M_{\kappa'} } } \sum_{l'} \Phi_{ \mu \mu' }(0\kappa, l'\kappa') e^{i \mathbf{q} \cdot \mathbf{r}(l')}.
    ```

## Action on displacements

We define left group action for positions {math}`\mathbf{r}(l\kappa) = \mathbf{r}(l) + \mathbf{r}(0\kappa)` by {math}`g = (\mathbf{R}_{g}, \mathbf{\tau}_{g}) \in \mathcal{G}` as
```{math}
  g \mathbf{r}(l\kappa)
    := \mathbf{R}_{g} \mathbf{r}(l\kappa) + \mathbf{\tau}_{g}.
```
We denote that site {math}`\mathbf{r}(0, \kappa)` is transformed to {math}`\mathbf{r}(0, g \kappa ) + \mathbf{h}_{g}(\kappa)` by symmetry operation {math}`g`.
Then,
```{math}
    g \mathbf{r}(\mathbf{l}, \kappa) &= \mathbf{R}_{g} \mathbf{r}(l) + \mathbf{r}(0, g\kappa) + \mathbf{h}_{g}(\kappa) \\
```
% \mathbf{h}_{g^{-1}}(\kappa) &= - \mathbf{p}_{g}^{-1} \mathbf{h}_{g}(g^{-1}\kappa)

We define left group action for displacement at {math}`\mathbf{r}(l\kappa)` as [^displacement_action]
```{math}
  g u_{\mu}(\mathbf{r}(l\kappa))
    := \sum_{\nu} [\mathbf{R}_{g}]_{\mu\nu} u_{\nu}(g \mathbf{r}(l\kappa)).
```

[^displacement_action]: This definition actually satisfies the condition of left group action
    ```{math}
        \left[ g \left( g' \mathbf{u}(\mathbf{r}(l\kappa)) \right) \right]_{\mu}
        &= \left[ g \left\{ \sum_{\nu} R_{g', \mu'\nu} u_{\nu}( g'\mathbf{r}(l\kappa) ) \right\}_{\mu'} \right]_{\mu} \\
        &= \sum_{ \mu'\nu } R_{g, \mu\mu'} R_{g', \mu'\nu} u_{\nu}( gg'\mathbf{r}(l\kappa) ) \\
        &= \sum_{ \nu } R_{gg', \mu\nu} u_{\nu}( gg'\mathbf{r}(l\kappa) ) \\
        &= \left[ (gg') \mathbf{u}(\mathbf{r}(l\kappa)) \right]_{\mu}.
    ```

Consider Fourier transformation of {math}`\mathbf{u}(\mathbf{r}(l\kappa))`
```{math}
  \mathbf{u}(\kappa; \mathbf{q})
    &:= \frac{1}{\sqrt{L^{3}}} \sum_{l} \mathbf{u}(\mathbf{r}(l\kappa)) e^{ i \mathbf{q} \cdot \mathbf{r}(l) } \\
  \mathbf{u}(\mathbf{r}(l\kappa))
    &= \frac{1}{\sqrt{L^{3}}} \sum_{\mathbf{q}} \mathbf{u}(\kappa; \mathbf{q}) e^{ -i \mathbf{q} \cdot \mathbf{r}(l) } \\
  g u_{\mu}(\kappa; \mathbf{q})
    &= \frac{1}{\sqrt{L^{3}}} \sum_{l} g u_{\mu}(\mathbf{r}(l\kappa)) e^{i \mathbf{q}\cdot \mathbf{r}(l)} \\
    &= \frac{1}{\sqrt{L^{3}}} \sum_{l}\sum_{\nu} R_{g,\mu\nu} u_{\nu}\left( \mathbf{R}_{g}\mathbf{r}(l) + \mathbf{r}(0, g\kappa) + \mathbf{h}_{g}(\kappa) \right) e^{i \mathbf{q}\cdot \mathbf{r}(l)} \\
    &= \frac{1}{\sqrt{L^{3}}} \sum_{l'}\sum_{\nu} R_{g,\mu\nu} u_{\nu}\left( \mathbf{r}(l', g\kappa) \right) e^{i \mathbf{q}\cdot \mathbf{R}_{g}^{-1}(\mathbf{r}(l') - \mathbf{h}_{\kappa} ) } \\
    &= \frac{1}{\sqrt{L^{3}}} \sum_{l'}\sum_{\nu} R_{g,\mu\nu} u_{\nu}\left( \mathbf{r}(l', g\kappa) \right) e^{i \mathbf{R}_{g}\mathbf{q}\cdot (\mathbf{r}(l') - \mathbf{h}_{\kappa} ) }
        \quad (\because \mathbf{R}_{g} \in O(3)) \\
    &= \sum_{\kappa'\mu'} u_{\mu'}(\kappa'; \mathbf{R}_{g} \mathbf{q} ) \Gamma_{\kappa'\mu'; \kappa\mu}^{\mathbf{q}}(g),
```
where
```{math}
:label: dynamical_matrix_rep
  \Gamma_{\kappa'\mu'; \kappa\mu}^{\mathbf{q}}(g)
    := \exp \left( -i \mathbf{R}_{g} \mathbf{q} \cdot \mathbf{h}_{g}(\kappa) \right) [\mathbf{R}_{g}]_{\mu'\mu} \delta_{ g\kappa, \kappa' }.
```
Equation {eq}`dynamical_matrix_rep` is essentially the same with Eq. (2.37) of {cite}`RevModPhys.40.1`.

We write matrix representation {math}`[\mathbf{\Gamma}^{\mathbf{q}}(g)]_{ \kappa'\mu'; \kappa\mu } := \Gamma_{\kappa'\mu'; \kappa\mu}^{\mathbf{q}}(g)`.
Then,
```{math}
    \left[ \mathbf{\Gamma}^{ \mathbf{q}}(gg') \right]_{ \kappa'\mu', \kappa\mu }
    &= \delta_{gg'\kappa, \kappa'} R_{gg', \mu'\mu} \exp \left( -i \mathbf{R}_{gg'}\mathbf{q} \cdot \mathbf{h}_{gg'}(\kappa) \right) \\
    &= \sum_{ \kappa''\nu }
        \delta_{g\kappa'', \kappa'} \delta_{g'\kappa, \kappa''}
        R_{g, \mu'\nu} R_{g', \nu\mu}
        \exp \left( -i \mathbf{R}_{g'}\mathbf{q} \cdot \mathbf{h}_{g'}(\kappa) \right)
        \exp \left( -i \mathbf{R}_{g}\mathbf{R}_{g'}\mathbf{q} \cdot \mathbf{h}_{g}(g\kappa') \right) \\
    &= \left[ \mathbf{\Gamma}^{ \mathbf{R}_{g'} \mathbf{q}}(g) \mathbf{\Gamma}^{ \mathbf{q}}(g') \right]_{ \kappa'\mu', \kappa\mu } \\
  \mathbf{\Gamma}^{\mathbf{q}}(g)^{\dagger} \mathbf{\Gamma}^{ \mathbf{q}}(g)
    &= \mathbf{1} \quad \mbox{(Unitary)}
```

Fourier transformation of force constants
```{math}
  \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q})
    &:= \sum_{l'} \Phi_{\mu\mu'}(0\kappa; l'\kappa') e^{ i \mathbf{q} \cdot \mathbf{r}(l') } \\
  \sum_{ l l' } \Phi_{ \mu \mu' }(l\kappa, l'\kappa') u_{\mu}(l\kappa) u_{\mu'}(l'\kappa')
    &= \sum_{\mathbf{q}} \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) u_{\mu}(\kappa; \mathbf{q}) u_{\mu'}(\kappa'; -\mathbf{q}) \\
```

The condition that potential energy is invariant under symmetry operations is rewritten as [^fourier_force_constant]
```{math}
  \mathbf{\Phi}(\mathbf{R}_{g} \mathbf{q})
    = \mathbf{\Gamma}^{\mathbf{q}}(g) \mathbf{\Phi}(\mathbf{q}) \mathbf{\Gamma}^{\mathbf{q}}(g)^{\dagger}.
```

[^fourier_force_constant]: The derivation is as follows:
    ```{math}
    \sum_{ l\kappa\mu l'\kappa'\mu' } \Phi_{ \mu \mu' }(l\kappa, l'\kappa') u_{\mu}(l\kappa) u_{\mu'}(l'\kappa')
    &= \sum_{ l\kappa\mu l'\kappa'\mu' } \Phi_{ \mu \mu' }(l\kappa, l'\kappa') gu_{\mu}(l\kappa) gu_{\mu'}(l'\kappa') \\
    &= \dots = \sum_{  } \left[ \mathbf{\Gamma}^{\mathbf{q}}(g) \mathbf{\Phi}(\mathbf{q}) \mathbf{\Gamma}^{\mathbf{q}}(g)^{\dagger} \right]_{\kappa\mu, \kappa'\mu'} u_{\mu}(\kappa; \mathbf{R}_{g}\mathbf{q}) u_{\mu'}(\kappa'; -\mathbf{R}_{g}\mathbf{q}).
    ```

## Diagonalizing {math}`\mathbf{\Phi}(\mathbf{q})` with modified eigenvectors

For {math}`h, h' \in \mathcal{G}^{\mathbf{q}}`,
```{math}
  \mathbf{\Gamma}^{ \mathbf{q}}(h) \mathbf{\Gamma}^{ \mathbf{q}}(h')
    &= \mathbf{\Gamma}^{ \mathbf{q}}(hh') \\
  \mathbf{\Phi}(\mathbf{q})
    &= \mathbf{\Gamma}^{\mathbf{q}}(h) \mathbf{\Phi}(\mathbf{q}) \mathbf{\Gamma}^{\mathbf{q}}(h)^{\dagger}
    \quad (\forall h \in \mathcal{G}^{\mathbf{q}}).
```
We can reuse eigenvectors of dynamical matrix for decomposing {math}`\mathbf{\Gamma}^{ \mathbf{q}}` into irreps,
```{math}
  D_{\mu\mu'}(\kappa\kappa'; \mathbf{q})
    &= \frac{1}{\sqrt{ M_{\kappa}M_{\kappa'} }} e^{ i \mathbf{q} \cdot \left( \mathbf{r}(0\kappa') - \mathbf{r}(0\kappa) \right) } \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) \\
  f_{\mu}(\kappa; \mathbf{q}\nu)
    &:= \frac{1}{\sqrt{M_{\kappa}}} e^{ i\mathbf{q} \cdot \mathbf{r}(0\kappa) } e_{\mu}(\kappa; \mathbf{q}\nu) \\
  u_{\mu}(l\kappa)
    &= \frac{1}{\sqrt{L^{3}}} \sum_{ q } Q_{q} f_{\mu}(\kappa; q) e^{ i \mathbf{q} \cdot \mathbf{r}(l) } \\
  \sum_{ \kappa \mu }\sum_{ \kappa' \mu' } f_{\mu}(\kappa; \mathbf{q}\nu)^{\ast} \Phi_{\mu\mu'}(\kappa\kappa'; \mathbf{q}) f_{\mu'}(\kappa'; \mathbf{q}\nu')
    &= \omega_{\mathbf{q}\nu}^{2} \delta_{\nu\nu'}
```

*modified eigenvectors*, {math}`\mathbf{f}(\mathbf{q}\nu) = \{f_{\mu}(\kappa; \mathbf{q}\nu) \}_{\kappa\mu}`, with the same eigenvalues form irreps,
```{math}
    \mathbf{F}_{\mathbf{q} \omega}
        &= \left( \mathbf{f}(\mathbf{q}\nu) \mid \omega_{\mathbf{q}\nu}^{2} = \omega^{2} \right) \\
    \mathbf{\Gamma}^{\mathbf{q}}_{\omega}
        &:= \mathbf{F}_{\mathbf{q} \omega}^{\dagger} \mathbf{\Gamma}^{\mathbf{q}} \mathbf{F}_{\mathbf{q} \omega}
```

## Small represenation of {math}`\mathcal{G}^{\mathbf{q}}`

We can introduce projective representation {math}`\mathbf{\gamma}^{ \mathbf{q}}(h)`,
```{math}
  \mathbf{\Gamma}^{ \mathbf{q}}(h)
    &=: e^{ -i \mathbf{q} \cdot \mathbf{v}_{h} } \mathbf{\gamma}^{ \mathbf{q}}(h) \\
  \mathbf{\gamma}^{ \mathbf{q}}(h) \mathbf{\gamma}^{ \mathbf{q}}(h')
    &= e^{ -i \mathbf{q} \cdot ( \mathbf{R}_{h} \mathbf{v}_{h'} - \mathbf{v}_{h'} ) } \mathbf{\gamma}^{ \mathbf{q}}(hh') \\
  \mathbf{\gamma}^{ \mathbf{q}}((E, \mathbf{t}))
    &= \mathbf{1}
```

Thus, we only need to consider projective representation {math}`\mathbf{\gamma}^{ \mathbf{q}}` for little co-group {math}`\overline{\mathcal{G}}^{\mathbf{q}} \simeq \mathcal{G}^{\mathbf{q}} / \mathcal{T}`.
We obtain projective irreps
```{math}
    \mathbf{\gamma}^{\mathbf{q}\omega}(\mathbf{R}_{h})
        := e^{ -i \mathbf{q} \cdot \mathbf{v}_{h} } \mathbf{F}_{\mathbf{q} \omega}^{\dagger} \mathbf{\gamma}^{\mathbf{q}}(\mathbf{R}_{h}) \mathbf{F}_{\mathbf{q} \omega}
        \quad (h \in \mathcal{G}^{\mathbf{q}})
```

## References

```{bibliography}
:filter: docname in docnames
```
