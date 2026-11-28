# Choose a representation theory topic

Find the explanation for the vector space and symmetry operations in your calculation.

!!! info "Before you start"

    You need linear algebra and crystallographic symmetry, plus the
    [spgrep conventions](../../conventions.md). These pages develop the theory;
    the [core API reference](../../api/api_core.md) lists callable functions.

Read [space-group irreps](spacegroup_irreps.md) first for translation phases,
little groups, and projective representations. Then choose the relevant extension:

- [Physically irreducible representations](reality.md) explains real vector spaces
  and the role of the pair of wave vectors $\pm\mathbf{k}$.
- [Spinor representations](spinor.md) explains the spinor factor system and rotation conventions.
- [Co-representations](corep.md) explains anti-linear symmetry operations and magnetic groups.

Once the representation type is clear, read [enumeration algorithms](../algorithm/algorithm.md)
or [symmetry-adapted bases](../projection.md) for how spgrep calculates and uses it.
