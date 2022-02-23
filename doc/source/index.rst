###################
ComCH documentation
###################

ComCH is a Python 3 package for the study of commutativity up to coherent homotopies.


Motivation
----------

Commutativity up to coherent homotopies is a concept originated in algebraic topology. It has found modern uses in topological data analysis, motion planning, condensed matter physics and several other areas. An important challenge for the application of the mathematical ideas surrounding this concept, which are often defined non-constructively, is to describe them in effective terms suitable for concrete computations. This package, a specialized computer algebra system, serves to bridge this gap between theoretical concepts and concrete applications.


Mathematical overview
---------------------

After the pioneering work of Steenrod, Adem, Serre, Cartan, Araki-Kudo, Dyer-Lashof, Stashef, Boardman-Vogt, May, and many others, today there is a rich theory of commutativity up to coherent homotopies whose modern framework is provided by operads and PROPs, and where :math:`E_n`-operads play a central role parameterizing the different levels of homotopical commutativity. In this package, we focus on the category of chain complexes, and consider two models of the :math:`E_\infty` equipped with filtrations by :math:`E_n`-operads.
These are respectively due to McClure-Smith and Berger-Fresse and are known as the surjection and Barratt-Eccles operads.


Installation
------------

This package is written in Python 3 and has no dependencies. It can be installed from a terminal simply entering:

:code:`python3 -m pip install comch`


.. toctree::
   :caption: API reference
   :maxdepth: 1
   
   modules/index


.. toctree::
   :caption: Jupyter notebooks
   :maxdepth: 1
   
   notebooks/index


Home
----

https://github.com/ammedmar/comch


References
----------

[McS]: J. McClure, and J. Smith. "Multivariable cochain operations and little n-cubes." Journal of the American Mathematical Society 16.3 (2003): 681-704.

[BF]: C. Berger, and B. Fresse. "Combinatorial operad actions on cochains." Mathematical Proceedings of the Cambridge Philosophical Society. Vol. 137. No. 1. Cambridge University Press, 2004.

[KMM]: R. Kaufmann, and A. Medina-Mardones. "Cochain level May-Steenrod operations." Forum Mathematicum 33 (2021), no. 6, 1507--1526.
