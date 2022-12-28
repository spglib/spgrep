---
title: 'spgrep: On-the-fly generator of space-group irreducible representations'
tags:
  - Python
  - computational materials science
  - crystallography
  - group theory
  - representation theory
  - irreducible representation
authors:
  - name: Kohei Shinohara
    orcid: 0000-0002-5907-2549
    corresponding: true
    affiliation: 1
  - name: Atsushi Togo
    orcid: 0000-0001-8393-9766
    affiliation: "2, 3"
  - name: Isao Tanaka
    orcid: 0000-0002-4616-118X
    affiliation: "1, 3, 4"
affiliations:
  - name: Department of Materials Science and Engineering, Kyoto University, Sakyo, Kyoto, Japan
    index: 1
  - name: Research and Services Division of Materials Data and Integrated System, National Institute for Materials Science, Tsukuba, Ibaraki, Japan
    index: 2
  - name: Center for Elements Strategy Initiative for Structural Materials, Kyoto University, Sakyo, Kyoto, Japan
    index: 3
  - name: Nanostructures Research Laboratory, Japan Fine Ceramics Center, Nagoya, Japan
    index: 4
date: 19 December 2022
bibliography: paper.bib
---

# Summary
<!-- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

The group theory and representation theory provide a formal and helpful way to exploit the symmetry of systems in condensed-matter physics and materials science [@Inui1996-et; @el2008symmetry; @Bradley2009-ze; @dresselhaus2010group].
When we consider the microscopic structure of a crystal, its symmetry is classified by space groups [@ITA2016].
Irreducible representations (irreps) of space groups serve as fundamental building blocks for classifying physical states and simplifying numerical calculations for crystals.
Although irreps of space groups were tabulated in seminal works [@faddeyev1964; @kovalev1965irreducible; @miller1967tables; @Zak1969; @Bradley2009-ze; @cracknell1979kronecker; @altmann1994point], it is tedious and error-prone to look up these tables.
`spgrep` is a Python package to enumerate irreps of space groups from given crystal structures or symmetry operations without the need for consulting databases of irreps.
`spgrep` computes various kinds of representations: (1) linear irreps, (2) physically irreducible irreps [@PhysRevB.43.11010], (3) projective irreps for spinor, and (4) projective irreducible co-representations for spinor.
`spgrep` can also construct symmetry-adapted bases of the obtained irreps.

<!--
- Enumerate irreps of crystallographic point groups as well
-->

# Statement of need
<!-- A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work. -->

<!-- Related packages tabulating irreps in BC -->
There are several packages or services to provide irreps of space groups originally tabulated in @Bradley2009-ze, including `ISOTROPY` Software Suite [@Stokes:pc5025; @ISOTROPY], Bilbao Crystallographic Server [@aroyo2011crystallography; @AroyoPerezMatoCapillasKroumovaIvantchevMadariagaKirovWondratschek+2006+15+27; @Aroyo:xo5013; @Elcoro:ks5574], `SpaceGroupIrep` [@LIU2021107993], and `IrRep` [@IRAOLA2022108226].
While these can be accessed from programs, they still require the user to follow the convention of settings for space groups.
In contract, `spgrep` only requires minimal knowledge for space groups because it computes irreps on the fly.
Additionally, `spgrep` can be used in conjunction with these tabulation-based packages such as when we need to assign historical labels for irreps.

<!-- Related packages computing characters from eigenvectors -->
For density functional theory calculations, there are a few packages to compute irreps from Bloch wave functions, including `Irvsp` [@GAO2021107760] and `qeirreps` [@MATSUGATANI2021107948].
Although these packages do not rely on tabulations of irreps, an arbitrary irrep will be obtained within unitary equivalence for multi-dimensional irreps, which undermines a deterministic symmetry-adapted basis.
On the other hand, `spgrep` provides unique irreps for given space groups through the implementation of a deterministic algorithm in @Neto:a09740.

<!--
Intertwiner on the fly [@THOMAS201776; @doi:10.1137/090779966]
-->

<!-- Minimal dependency and permissive license -->
It is advantageous for a package for irreps to easily integrate with other domain-specific packages in condensed-matter physics and materials science.
With this in mind, `spgrep` has been implemented with minimal dependencies: `numpy` [@harris2020array] for array programming and `spglib` [@spglibv1; @spglibv2] for crystal symmetry search.
Also, `spgrep` is distributed under the permissive BSD 3-clause license.

<!-- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it. -->
In particular, `spgrep` is well-suited for implementation of group-theoretic methods with automated calculations.
We plan to implement a group-theoretic method for analyzing continuous phase transition, called isotropy subgroup [@doi:10.1142/0751], on top of `spgrep` and apply it to automated phase transition path search [@spgrep_application].

# Acknowledgements
<!-- Acknowledgement of any financial support. -->

This work was supported by a Grant-in-Aid for JSPS Research Fellows (Grant Number 21J10712) from the Japan Society for the Promotion of Science (JSPS).

# References
<!-- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline. -->
