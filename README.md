[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ammedmar/comch/master?filepath=notebooks)

# ComCH

ComCH is a Python 3 package for the study of commutativity
up-to-coherent-homotopies.

## Motivation

Commutativity up-to-coherent-homotopies is a concept originated in
algebraic topology. It has found modern uses in topological data
analysis, motion planning, condensed matter physics and several other
areas. An important challenge for the application of the mathematical
ideas sourraounding this concept, which are often defined
non-constructively, is to describe them in effective terms suitable for
concrete computations. This package serves to bridge the gap between
theoretical ideas and concrete applications, by implementing the
structures needed to interact effectively with commutativity
up-to-coherent-homotopies.

## Mathematical overview

Following the pioneering work of Steenrod, Cartan, Adem, Stashef,
Bordman-Vogt, May, Dyer-Lashof and others, today we understand the
correct framework for the study of commutativity
up-to-coherent-homotopies as the one provided by operads and PROPs. In
particular, *E*<sub>*n*</sub>-operads play a central role parameterizing
the different levels of homotopical commutativity. In this package, we
focus on the category of chain complexes, and consider two models for
the *E*<sub>∞</sub>-operad which are equipped with a filtrations by
*E*<sub>*n*</sub>-operads. These models are respectively due to
McClure-Smith [McS] and Berger-Fresse [BF] and are known as the
surjection and Barratt-Eccles operads.

## Installation

This package is written in Python 3 and has no dependencies. It can be
installed from a terminal simply entering:

`python3 -m pip install comch`

## References

[McS]: J. McClure, and J. Smith. "Multivariable cochain operations and
little n-cubes." Journal of the American Mathematical Society 16.3
(2003): 681-704.

[BF]: C. Berger, and B. Fresse. "Combinatorial operad actions on
cochains." Mathematical Proceedings of the Cambridge Philosophical
Society. Vol. 137. No. 1. Cambridge University Press, 2004.
