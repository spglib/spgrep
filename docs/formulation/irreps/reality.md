(physically_irreps)=
# Physically irreducible representation (PIR)

Let $(\Gamma, \mathrm{Span}_{\mathbb{C}} \{ \mathbf{v}_{i} \}_{i=1}^{d} )$ be a unitary (projective) irrep of group $G$.
Then, its conjugate (projective) representation is $(\Gamma^{\ast}, \mathrm{Span}_{\mathbb{C}} \{ \mathbf{v}_{i}^{\ast} \}_{i=1}^{d} )$.
There are three cases for $\Gamma$ and $\Gamma^{\ast}$:
  1. $\Gamma$ is real: $\Gamma$ and $\Gamma^{\ast}$ are equivalent and can be taken as real matrices.
  2. $\Gamma$ is pseudo-real: $\Gamma$ and $\Gamma^{\ast}$ are equivalent but can be taken as real matrices.
  3. $\Gamma$ is not equivalent to $\Gamma^{\ast}$.

We sometimes need to restrict irrep under a vector space over $\mathbb{R}$ (instead of $\mathbb{C}$), which is called physically irreducible representation (PIR) {cite}`PhysRevB.43.11010`.

## Frobenius-Schur indicator

### Finite group

These cases are classified with Frobenius-Schur indicator:
$$
  \frac{1}{|G|} \sum_{ g \in G } \chi(g^{2})
  = \begin{cases}
    1 & \mbox{($\Gamma$ is real)} \\
    -1 & \mbox{($\Gamma$ is pseudo-real)} \\
    0 & \mbox{($\Gamma$ is not equivalent to $\Gamma^{\ast}$)} \\
  \end{cases}.
$$

### Space-group representations

Let $\Gamma^{(\mathbf{k}, 0)}$ and $\Gamma^{(\mathbf{k}, 1)}$ be representations of space group $\mathcal{G}$ with $\mathbf{k}$ vector.
Let $\mathcal{T}$ be a translational subgroup of $\mathcal{G}$.
Consider coset representatives of little group of $\mathbf{k}$ over $\mathcal{T}$:
$$
  \mathcal{G}^{\mathbf{k}} = \coprod_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } (\mathbf{S}_{i}, \mathbf{w}_{i}) \mathcal{T}.
$$

Then sum over translations in $\mathcal{G}^{\mathbf{k}}$:
$$
  \frac{1}{|\mathcal{G}^{\mathbf{k}}|} \sum_{ g \in \mathcal{G}^{\mathbf{k}} } \chi^{\mathbf{k}}(g^{2})
  &= \frac{1}{N} \sum_{ \mathbf{t} } \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
      \chi^{\mathbf{k}}\left( (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) \right) \\
  &= \frac{1}{N} \sum_{ \mathbf{t} } \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
      \chi^{\mathbf{k}}\left( (\mathbf{E}, \mathbf{t} + \mathbf{S}_{i}\mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i})^{2} \right) \\
  &= \left(
        \frac{1}{N} \sum_{ \mathbf{t} } e^{-i \mathbf{k} \cdot (\mathbf{t} + \mathbf{S}_{i}\mathbf{t}) }
     \right)
     \left(
        \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } \chi^{\mathbf{k}}\left( (\mathbf{S}_{i}, \mathbf{w}_{i})^{2} \right)
     \right) \\
  &= \mathbb{I}\left[ 2\mathbf{k} \equiv \mathbf{0} \right] \cdot
      \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} } \chi^{\mathbf{k}}\left( (\mathbf{S}_{i}, \mathbf{w}_{i})^{2} \right),
$$
where $\mathbb{I}[C]$ takes one if the condition $C$ is true and takes zero otherwise [^derivation-2k].

[^derivation-2k]: Here we use $\mathbf{S}_{i}^{\top} \in \overline{\mathcal{G}}^{\mathbf{k}}$ as
  $$
    \frac{1}{N} \sum_{ \mathbf{t} } e^{-i \mathbf{k} \cdot (\mathbf{t} + \mathbf{S}_{i}\mathbf{t}) }
      &= \frac{1}{N} \sum_{ \mathbf{t} } e^{-i (\mathbf{E} + \mathbf{S}_{i})^{\top} \mathbf{k} \cdot \mathbf{t} } \\
      &= \mathbb{I} \left[ (\mathbf{E} + \mathbf{S}_{i})^{\top} \mathbf{k} \equiv \mathbf{0} \right] \\
      &= \mathbb{I} \left[ \mathbf{k} + \mathbf{S}_{i}^{\top} \mathbf{k} \equiv \mathbf{0} \right] \\
      &= \mathbb{I} \left[ 2\mathbf{k} \equiv \mathbf{0} \right].
  $$

## PIR of finite group

### (1) $\Gamma$ is real

The following construction is based on Ref. {cite}`Inui1996-et`.

In this case, since $\Gamma$ and $\Gamma^{\ast}$ are equivalent, there exists a symmetric unitary matrix with
$$
  \mathbf{\Gamma}(g) \mathbf{U} &= \mathbf{U} \mathbf{\Gamma}(g)^{\ast} \\
  \mathbf{U}^{\dagger} \mathbf{U} &= \mathbf{1} \\
  \mathbf{U}^{\top} &= \mathbf{U}.
$$
We impose $\det \mathbf{U} = 1$ additionally to specify the intertwiner uniquely.
The intertwiner can be numerically computed as shown in {ref}`here <intertwiner>`.

The symmetric unitary matrix $\mathbf{U}$ can be diagonalized with real orthogonal matrix $\mathbf{S}$ as $\mathbf{U} = \mathbf{S}^{-1} \mathbf{\Omega} \mathbf{S}$ [^sym_unitary].
For this symmetric unitary matrix $\mathbf{U}$, we can choose its square root with
$$
  \mathbf{T} &:= \mathbf{S}^{-1} \mathbf{\Omega}^{1/2} \mathbf{S} \\
  \mathbf{T}^{2} &= \mathbf{U} \\
  \mathbf{T}^{\dagger} \mathbf{T} &= \mathbf{1} \\
  \mathbf{T}^{\top} &= \mathbf{T}.
$$

A transformed representation $\mathbf{\Gamma}'(g) := \mathbf{T}^{-1}\mathbf{\Gamma}(g)\mathbf{T}$ is real because
$$
  \mathbf{\Gamma}'(g)
  &= \mathbf{T}^{-1} \mathbf{U} \mathbf{\Gamma}(g)^{\ast} \mathbf{U}^{-1} \mathbf{T} \\
  &= \mathbf{T} \mathbf{\Gamma}(g)^{\ast} \mathbf{T}^{-1} \\
  &= \mathbf{\Gamma}'(g)^{\ast}.
$$

[^sym_unitary]: For symmetric unitary matrix $\mathbf{U}$, if $\mathbf{v}$ is eigenvector of $\mathbf{U}$, its conjugacy $\mathbf{v}'$ is also eigenvector.
    Thus, we can take real part or imaginary part of $\mathbf{v}$ as a new basis vector.

### (2, 3) $\Gamma$ is pseudo-real or not equivalent to $\Gamma^{\ast}$

In these cases, we can transform conjugated basis pair to real vectors by unitary matrix:

$$
  (\mathbf{v}_{1}, \cdots, \mathbf{v}_{d}, \mathbf{v}_{1}^{\ast}, \cdots, \mathbf{v}_{d}^{\ast}) \mathbf{U}
    &= \sqrt{2} (\mathrm{Re}\, \mathbf{v}_{1}, \cdots, \mathrm{Re}\, \mathbf{v}_{d}, \mathrm{Im}\, \mathbf{v}_{1}, \cdots, \mathrm{Im}\, \mathbf{v}_{d}) \\
  \mathbf{U} &:= \frac{1}{\sqrt{2}}\begin{pmatrix}
    \mathbf{1}_{d} & -i \mathbf{1}_{d} \\
    \mathbf{1}_{d} & i \mathbf{1}_{d} \\
  \end{pmatrix} \quad (\mathrm{Unitary}) \\
  \mathbf{U}^{-1}
  \begin{pmatrix}
    \mathbf{D}(g) & \\
    & \mathbf{D}(g)^{\ast}
  \end{pmatrix}
  \mathbf{U}
  &= \begin{pmatrix}
    \mathrm{Re}\, \mathbf{D}(g) & \mathrm{Im}\, \mathbf{D}(g) \\
    -\mathrm{Im}\, \mathbf{D}(g) & \mathrm{Re}\, \mathbf{D}(g) \\
  \end{pmatrix}
  \quad (g \in G)
$$

## PIR of space group $\mathcal{G}$

To construct PIR of space group from small representation $\Gamma^{\mathbf{k}\alpha}$ at $\mathbf{k}$ {cite}`Stokes:pc5025`, we apply the finite-group recipe above to the little group $\mathcal{G}^{\mathbf{k}}$.
So that the construction remains valid when $2\mathbf{k} \not\equiv \mathbf{0}$, spgrep (as of v0.6.0) enlarges the carrier from $\mathcal{G}^{\mathbf{k}}$ to the *little group $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}$ of $\pm\mathbf{k}$* following {cite}`Stokes:pc5025`; this matches the presentation used by the Bilbao Crystallographic Server.

(pir_kbark)=
### Little group $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}$ of $\pm\mathbf{k}$

Define the *little co-group of $\pm\mathbf{k}$* as the set of point-group operations that stabilize $\{ \mathbf{k}, -\mathbf{k} \}$ as a set:

$$
  \overline{\mathcal{G}}^{\mathbf{k}\bar{\mathbf{k}}}
  := \left\{ \mathbf{S} \in \overline{\mathcal{G}} \;\middle|\; \mathbf{S}^{\top} \mathbf{k} \equiv \pm \mathbf{k} \pmod{\mathbf{K}} \right\}.
$$

By construction, any $\mathbf{S}$ stabilizing $\mathbf{k}$ also stabilizes $\{ \mathbf{k}, -\mathbf{k} \}$, so $\overline{\mathcal{G}}^{\mathbf{k}} \subseteq \overline{\mathcal{G}}^{\mathbf{k}\bar{\mathbf{k}}} \subseteq \overline{\mathcal{G}}$.
The *little group of $\pm\mathbf{k}$* is its preimage under the point-group projection,

$$
  \mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}
  := \left\{ g = (\mathbf{S}_{g}, \mathbf{w}_{g}) \in \mathcal{G} \;\middle|\; \mathbf{S}_{g} \in \overline{\mathcal{G}}^{\mathbf{k}\bar{\mathbf{k}}} \right\}
  = \coprod_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}\bar{\mathbf{k}}} \} } (\mathbf{S}_{i}, \mathbf{w}_{i}) \mathcal{T},
$$

with $\mathcal{G}^{\mathbf{k}} \subseteq \mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}} \subseteq \mathcal{G}$.
When $2\mathbf{k} \equiv \mathbf{0}$ we have $\mathbf{k} \equiv -\mathbf{k}$, so $\overline{\mathcal{G}}^{\mathbf{k}\bar{\mathbf{k}}} = \overline{\mathcal{G}}^{\mathbf{k}}$ and $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}} = \mathcal{G}^{\mathbf{k}}$ automatically, so $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}$ gives a uniform domain regardless of whether $2\mathbf{k} \equiv \mathbf{0}$.

### Real form on $\mathcal{G}^{\mathbf{k}}$

**Case (1): $\Gamma^{\mathbf{k}\alpha}$ is real.**
Averaging the finite-group intertwiner of {ref}`(1) <physically_irreps>` over translations, and using $2\mathbf{k} \equiv \mathbf{0}$ (which holds whenever $\Gamma^{\mathbf{k}\alpha}$ is real) to collapse the translation sum via $2\mathbf{k} \cdot \mathbf{t} \in \mathbb{Z}$:
$$
    \mathbf{U}
    &:= \frac{1}{N} \sum_{ \mathbf{t} } \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) \right)
            \mathbf{B}
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{E}, \mathbf{t})(\mathbf{S}_{i}, \mathbf{w}_{i}) \right)^{\ast \dagger} \\
    &= \sum_{ \{ i \mid \mathbf{S}_{i} \in \overline{\mathcal{G}}^{\mathbf{k}} \} }
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)
            \mathbf{B}
            \mathbf{\Gamma}^{\mathbf{k}\alpha}\left( (\mathbf{S}_{i}, \mathbf{w}_{i}) \right)^{\ast \dagger}.
$$
Feed this $\mathbf{U}$ through the $\mathbf{T} = \mathbf{U}^{1/2}$ step of the finite-group (1) subsection to obtain the real form.

**Cases (2, 3): $\Gamma^{\mathbf{k}\alpha}$ is pseudo-real or not equivalent to $\Gamma^{\mathbf{k}\alpha \ast}$.**
Apply the block-diagonal identity from the finite-group (2, 3) subsection with $g \in \mathcal{G}^{\mathbf{k}}$.
For a translation $(\mathbf{E}, \mathbf{t}) \in \mathcal{G}^{\mathbf{k}}$ the block form reduces to a planar rotation:
$$
  \tilde{\mathbf{\Gamma}}^{\mathbf{k}\alpha}( (\mathbf{E}, \mathbf{t}) )
    = \begin{pmatrix}
      \cos (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} & -\sin (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} \\
      \sin (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} & \cos (\mathbf{k} \cdot \mathbf{t}) \mathbf{1}_{d} \\
    \end{pmatrix}.
$$

### Extension to $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}$

When $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}} \neq \mathcal{G}^{\mathbf{k}}$, fix any $a \in \mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}} \setminus \mathcal{G}^{\mathbf{k}}$ so that $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}} = \mathcal{G}^{\mathbf{k}} \sqcup a \mathcal{G}^{\mathbf{k}}$ (note $a^{2} \in \mathcal{G}^{\mathbf{k}}$).
Extending $\mathbf{\Gamma}^{\mathbf{k}\alpha}$ from $\mathcal{G}^{\mathbf{k}}$ to $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}$ is structurally identical to the co-representation construction in {ref}`corep <corep>`: take $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}$ in place of the magnetic group, $\mathcal{G}^{\mathbf{k}}$ in place of its unitary subgroup, $a$ as the antisymmetry coset representative, and let the role of complex conjugation be played by the anti-linear action $\mathbf{k} \to -\mathbf{k}$, so that the conjugated irrep is $\overline{\mathbf{\Gamma}}^{\mathbf{k}\alpha}(h) = \mathbf{\Gamma}^{\mathbf{k}\alpha}(a^{-1} h a)^{\ast}$ (trivial factor system).

The Frobenius-Schur indicator
$$
  \xi^{\alpha} := \frac{1}{|\mathcal{G}^{\mathbf{k}}|} \sum_{h \in \mathcal{G}^{\mathbf{k}}} \chi^{\mathbf{k}\alpha}((ah)^{2})
  \quad \in \{ -1, 0, 1 \}
$$
selects the block-matrix form of $\tilde{\mathbf{\Gamma}}^{\mathbf{k}\alpha}(a)$; see {ref}`corep <corep>` for the three cases and their derivation.
In every case $\tilde{\mathbf{\Gamma}}^{\mathbf{k}\alpha}(a h) = \tilde{\mathbf{\Gamma}}^{\mathbf{k}\alpha}(a) \tilde{\mathbf{\Gamma}}^{\mathbf{k}\alpha}(h)$ on the remaining coset $a \mathcal{G}^{\mathbf{k}}$.

### Realification of the corep on $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}$

Reuse the change of basis from the finite-group $(2, 3)$ subsection and write $\tilde{\mathbf{D}}(g) := \tilde{\mathbf{\Gamma}}^{\mathbf{k}\alpha}(g)$.
For $g \in \mathcal{G}^{\mathbf{k}}$ the corep acts $\mathbb{C}$-linearly, so on the doubled basis it is block-diagonal and the realification is the same formula as the finite-group $(2, 3)$ subsection:

$$
  \mathbf{U}^{-1}
  \begin{pmatrix}
    \tilde{\mathbf{D}}(g) & \\
    & \tilde{\mathbf{D}}(g)^{\ast} \\
  \end{pmatrix}
  \mathbf{U}
  = \begin{pmatrix}
    \mathrm{Re}\, \tilde{\mathbf{D}}(g) & \mathrm{Im}\, \tilde{\mathbf{D}}(g) \\
    -\mathrm{Im}\, \tilde{\mathbf{D}}(g) & \mathrm{Re}\, \tilde{\mathbf{D}}(g) \\
  \end{pmatrix}
  \quad (g \in \mathcal{G}^{\mathbf{k}}).
$$

For $g \in a\mathcal{G}^{\mathbf{k}}$ the corep acts $\mathbb{C}$-antilinearly as $\mathbf{v} \mapsto \tilde{\mathbf{D}}(g) \mathbf{v}^{\ast}$.
On the doubled basis $(\mathbf{v}_{1}, \cdots, \mathbf{v}_{d}, \mathbf{v}_{1}^{\ast}, \cdots, \mathbf{v}_{d}^{\ast})$ this antilinear action is the $\mathbb{C}$-linear block-off-diagonal matrix $\bigl(\begin{smallmatrix} \mathbf{0} & \tilde{\mathbf{D}}(g) \\ \tilde{\mathbf{D}}(g)^{\ast} & \mathbf{0} \end{smallmatrix}\bigr)$, since $g$ exchanges the $\mathbf{v}$ and $\mathbf{v}^{\ast}$ blocks of basis vectors while complex-conjugating coordinates.
Conjugating with $\mathbf{U}$ gives

$$
  \mathbf{U}^{-1}
  \begin{pmatrix}
    \mathbf{0} & \tilde{\mathbf{D}}(g) \\
    \tilde{\mathbf{D}}(g)^{\ast} & \mathbf{0} \\
  \end{pmatrix}
  \mathbf{U}
  = \begin{pmatrix}
    \mathrm{Re}\, \tilde{\mathbf{D}}(g) & -\mathrm{Im}\, \tilde{\mathbf{D}}(g) \\
    -\mathrm{Im}\, \tilde{\mathbf{D}}(g) & -\mathrm{Re}\, \tilde{\mathbf{D}}(g) \\
  \end{pmatrix}
  \quad (g \in a\mathcal{G}^{\mathbf{k}}).
$$

The corep multiplication law (with $\ast$ applied iff the left factor is antilinear) is reproduced by ordinary matrix multiplication of these $2d \times 2d$ real blocks, so the resulting matrices form a real linear representation of $\mathcal{G}^{\mathbf{k}\bar{\mathbf{k}}}$.

## References

```{bibliography}
:filter: docname in docnames
```
