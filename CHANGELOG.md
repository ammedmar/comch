# Changelog

All notable changes to ComCH are documented in this file.

## 0.2.0

### Added

- Chain-level representatives of Steenrod operations through the `Steenrod`
  namespace.
- Even-prime Cartan homotopies and simplicial Cartan coboundaries.
- Odd-prime Cartan relators, including the coproduct correction of
  Cantero-Morán and Medina-Mardones, and simplicial representatives in the
  Berger–Fresse sign convention.
- Even-prime Adem homotopies and complete Adem coboundaries, including the two
  non-symmetric correction terms.
- Direct action of Barratt–Eccles elements on simplicial chains.
- Guided Cartan and Adem notebooks with low-degree mathematical checks.
- Automated wheel and source-distribution publication through PyPI trusted
  publishing.

### Changed

- Renamed `may_steenrod_structure` to `may_steenrod_map` on `BarrattEccles`
  and `Surjection`.
- Renamed `steenrod_operation` to
  `steenrod_operation_representative` on those classes.
- Normalized torsion values and expanded support for the Berger–Fresse sign
  convention.

The former method names remain available in 0.2.0 and issue
`DeprecationWarning`; they are compatibility aliases and may be removed in a
future release.

### Known limitations

- The Cartan and Adem convenience functions currently evaluate on simplicial
  chains. The corresponding direct Barratt–Eccles action on cubical chains is
  not implemented.
- Odd-prime Cartan table reduction currently requires the Berger–Fresse sign
  convention.
