[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ammedmar/comch/master?urlpath=lab/tree/notebooks)

# ComCH

ComCH is a specialized computer algebra system for the study of commutativity
up to coherent homotopies. It provides effective constructions suitable for
concrete computations in algebraic topology.

See the [documentation](https://comch.readthedocs.io/en/latest/) for the API
reference and guided notebooks.

## Motivation

Commutativity up to coherent homotopies originates in algebraic topology and
has found modern uses in topological data analysis, motion planning, condensed
matter physics, and other areas. Many of the surrounding mathematical ideas
are defined non-constructively; ComCH helps bridge the gap between those
theoretical concepts and concrete applications.

## Mathematical overview

Following the pioneering work of Steenrod, Cartan, Adem, Serre, Araki-Kudo,
Dyer-Lashof, Stasheff, Boardman-Vogt, May, and others, the modern framework for
commutativity up to coherent homotopies is provided by operads and PROPs.
$E_n$-operads parameterize different levels of homotopical commutativity.
ComCH focuses on chain complexes and implements two models of the
$E_\infty$-operad equipped with filtrations by $E_n$-operads: the
McClure-Smith surjection operad [McS] and the Berger-Fresse Barratt-Eccles
operad [BF]. It also provides effective constructions of cochain-level
May-Steenrod operations [KMM].

## Installation

ComCH is written in pure Python and has no runtime dependencies:

```bash
python3 -m pip install comch
```

## Jupyter notebooks

To run the notebooks locally, install the notebook extra and register its
Python kernel:

```bash
python3 -m pip install "comch[notebooks]"
python3 -m ipykernel install --user --name comch --display-name "Python (comch)"
jupyter lab
```

Open an example notebook and select the `Python (comch)` kernel. The notebooks
can also be run directly in Binder using the badge at the top of this page.

## Development

Install the development extra from a repository checkout:

```bash
python3 -m pip install -e ".[dev]"
```

## References

[McS]: J. McClure, and J. Smith. "Multivariable cochain operations and little n-cubes." Journal of the American Mathematical Society 16.3 (2003): 681-704. [DOI](https://doi.org/10.1090/S0894-0347-03-00419-3)

[BF]: C. Berger, and B. Fresse. "Combinatorial operad actions on cochains." Mathematical Proceedings of the Cambridge Philosophical Society. Vol. 137. No. 1. Cambridge University Press, 2004. [DOI](https://doi.org/10.1017/S0305004103007138)

[KMM]: R. M. Kaufmann and A. M. Medina-Mardones. "Cochain level May-Steenrod operations." Forum Mathematicum 33 (2021), no. 6, 1507-1526. [DOI](https://doi.org/10.1515/forum-2020-0296)
