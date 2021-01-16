---

title: 'ComCH: A specialized Python algebra system for homotopical commutativity'

tags:

  - Python

  - computer algebra system

  - topology

  - homotopical algebra

  - operads

  - cohomology operations

authors:

  - name: Anibal M. Medina-Mardones

    orcid: 0000-0002-0905-262X

    affiliation: "1, 2" # (Multiple affiliations must be quoted)

affiliations:

 - name: Max Planck Institute for Mathematics, Bonn, Germany

   index: 1

 - name: University of Notre Dame, IN, USA

   index: 2

date: 11 November 2020

bibliography: paper.bib

---

# Summary

All the basic notions of number, from the integers to the complex, are equipped with a commutative product, and it was believed until Hamilton's introduction of the quaternions, that the product of any number system must be commutative. Hamilton's discovery allowed the consideration of other algebraic structures where commutativity is not assumed, and the conceptual shift that followed is only comparable to the effect non-euclidean geometries had in the study of shapes. Around a century later, after the development of the novel field of topology and homotopy, mathematicians returned to the question of commutativity and identified different levels laying in between the basic dichotomy. These intermediate structures correspond to coherent systems correcting homotopically the lack of strict commutativity, and constitute the focus of extensive research in mathematics.

# Statement of need

`ComCH` is a specialized computer algebra system for the study of commutativity up-to-coherent-homotopies, and it is built on `Python` with no dependencies. An important challenge for the application of commutativity up-to-coherent-homotopies and related notions, which are often defined non-constructively, is to describe them in effective terms suitable for concrete computations. This package serves to bridge this gap between theoretical concepts and concrete applications, by providing effective models for complexes parameterizing different levels of homotopical commutativity.

The theoretical concepts covered in this work have already found modern uses in topological data analysis, condensed matter physics, motion planning, and several other areas, see for example [@mm_persistence, @kapustin_thorngren, @g_lm_rm], with the notion of cohomology operations at the chain level [@mm, @kaufmann_mm, @brumfiel_mm_morgan] playing a central role.

# Mathematical overview

Following the pioneering work of Steenrod, Cartan, Adem, Stashef, Boardman-Vogt, May, Dyer-Lashof and others, today we understand the correct framework for the study of commutativity up-to-coherent-homotopies as the one provided by operads and PROPs. In particular, $E_n$-operads play a central role parameterizing the different levels of homotopical commutativity. In this package, we focus on the category of chain complexes, and consider two models for the $E_\infty$-operad which are equipped with filtrations by $E_n$-operads. These models are respectively due to McClure-Smith [@mcclure_smith] and Berger-Fresse [@berger_fresse] and are known as the surjection and Barratt-Eccles operads.

# Acknowledgements

We acknowledge contributions from Djian Post, Wojciech Reise and Michelle Smith, and the support of Kathryn Hess. Work partially supported by Innosuisse grant 32875.1 IP-ICT - 1.

# References