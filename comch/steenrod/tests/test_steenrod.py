import unittest
from itertools import combinations

from comch.simplicial import SimplicialElement
from comch.steenrod import Steenrod
from comch.surjection import Surjection, SurjectionElement


class TestSteenrod(unittest.TestCase):
    def test_chain_operations(self):
        x = Steenrod.chain_operations(2, -1, -2)
        y = SimplicialElement({
            ((0, 2, 3), (0, 1, 2)): 1,
            ((0, 1, 3), (1, 2, 3)): 1,
        }, torsion=2)
        self.assertEqual(x, y)

    def test_cartan_coboundaries(self):
        x = Steenrod.cartan_coboundaries(0, 3)
        y = SimplicialElement({
            ((0, 1), (1, 2), (1, 2), (2, 3)): 1,
        }, torsion=2)
        self.assertEqual(x, y)

        # the Barratt-Eccles action on cubical chains is not yet implemented
        with self.assertRaises(NotImplementedError):
            Steenrod.cartan_coboundaries(0, 3, shape='cube')

    def test_cartan_relation_on_cochains(self):
        r"""End-to-end check of the Cartan relation on the 4-simplex.

        For the mod-2 cocycles a = delta(dual of vertex 1) and
        b = delta(dual of vertex 3), checks on the top cell that

            (delta zeta_0)(a, b) =
            (a cup_0 b) cup_0 (a cup_0 b) + (a cup_0 a) cup_0 (b cup_0 b),

        where the cup products on the right are computed pointwise as
        cochain operations. Both sides equal 1, so the check is not vacuous.

        """
        n = 4
        cells = {d: list(combinations(range(n + 1), d + 1))
                 for d in range(n + 1)}

        def faces(cell):
            return [cell[:t] + cell[t + 1:] for t in range(len(cell))]

        def value(surjection, cochains, cell):
            chain = surjection(SimplicialElement({(cell,): 1}, torsion=2))
            return sum(v for k, v in chain.items()
                       if all(f in a for f, a in zip(k, cochains))) % 2

        cup_0 = SurjectionElement({(1, 2): 1}, torsion=2,
                                  convention='Berger-Fresse')
        zeta_0 = Surjection.cartan_homotopy(0)

        a = {(0, 1), (1, 2), (1, 3), (1, 4)}  # coboundary of dual of (1,)
        b = {(0, 3), (1, 3), (2, 3), (3, 4)}  # coboundary of dual of (3,)
        a_b = {c for c in cells[2] if value(cup_0, (a, b), c)}
        a_a = {c for c in cells[2] if value(cup_0, (a, a), c)}
        b_b = {c for c in cells[2] if value(cup_0, (b, b), c)}

        top = tuple(range(n + 1))
        rhs = (value(cup_0, (a_b, a_b), top)
               + value(cup_0, (a_a, b_b), top)) % 2
        lhs = sum(value(zeta_0, (a, a, b, b), face)
                  for face in faces(top)) % 2
        self.assertEqual(rhs, 1)
        self.assertEqual(lhs, rhs)

    def test_adem_coboundaries(self):
        x = Steenrod.adem_coboundaries(2, 1, 4)
        y = SimplicialElement({
            ((0, 1, 3), (1, 2, 3), (0, 3, 4), (0, 1, 3)): 1,
            ((0, 1, 3), (1, 2, 3), (0, 1, 4), (1, 3, 4)): 1,
            ((0, 2, 4), (2, 3, 4), (0, 1, 4), (1, 2, 4)): 1,
            ((0, 2, 4), (2, 3, 4), (0, 2, 4), (0, 1, 2)): 1,
            ((0, 1, 4), (0, 1, 4), (1, 2, 3), (1, 3, 4)): 1,
            ((0, 1, 4), (0, 1, 4), (2, 3, 4), (1, 2, 4)): 1,
            ((0, 1, 4), (0, 1, 4), (1, 3, 4), (1, 2, 3)): 1,
            ((0, 1, 4), (0, 1, 4), (1, 2, 4), (2, 3, 4)): 1,
            ((0, 3, 4), (0, 3, 4), (0, 2, 3), (0, 1, 2)): 1,
            ((0, 3, 4), (0, 3, 4), (0, 1, 3), (1, 2, 3)): 1,
            ((0, 3, 4), (0, 2, 3), (2, 3, 4), (0, 1, 2)): 1,
            ((0, 2, 3), (0, 1, 2), (0, 3, 4), (0, 2, 3)): 1,
            ((0, 3, 4), (1, 2, 3), (0, 1, 4), (1, 3, 4)): 1,
            ((0, 1, 2), (2, 3, 4), (0, 2, 3), (0, 1, 2)): 1,
            ((0, 1, 2), (2, 3, 4), (0, 3, 4), (0, 1, 3)): 1,
            ((0, 1, 2), (2, 3, 4), (0, 1, 3), (1, 2, 3)): 1,
            ((0, 1, 2), (2, 3, 4), (0, 1, 4), (1, 3, 4)): 1,
            ((0, 1, 2), (2, 3, 4), (1, 3, 4), (1, 2, 3)): 1,
            ((0, 1, 2), (2, 3, 4), (1, 2, 4), (2, 3, 4)): 1,
            ((0, 1, 2), (2, 3, 4), (0, 1, 4), (0, 1, 4)): 1,
            ((0, 1, 2), (2, 3, 4), (1, 2, 4), (1, 2, 4)): 1,
            ((0, 1, 2), (2, 3, 4), (2, 3, 4), (2, 3, 4)): 1,
        }, torsion=2)
        self.assertEqual(x, y)

        self.assertEqual(Steenrod.adem_coboundaries(1, 0, 3),
                         SimplicialElement(torsion=2))

        # the Barratt-Eccles action on cubical chains is not yet implemented
        with self.assertRaises(NotImplementedError):
            Steenrod.adem_coboundaries(1, 0, 1, shape='cube')

    def test_adem_relation_on_cochains(self):
        r"""End-to-end check of the Adem relation for (q, p) = (3, 1).

        For an explicit mod-2 cocycle alpha of degree 2 on the standard
        5-simplex, checks on a 4-cell where the left-hand side is nonzero
        that

            Sq^1 Sq^1 (alpha) = delta(A_{3,1}(alpha)),

        the (q, p) = (3, 1) instance of the cochain-level Adem relation of
        [BMM]. Here A_{3,1} is the action on alpha^{otimes 4} of the
        complete element returned by Surjection.adem_coboundary, including
        the homotopy component and both non-symmetric corrections.

        """

        def theta(i):
            key = tuple(1 if j % 2 == 0 else 2 for j in range(i + 2))
            return SurjectionElement({key: 1}, torsion=2,
                                     convention='Berger-Fresse')

        def value(surjection, cochains, cell):
            chain = surjection(SimplicialElement({(cell,): 1}, torsion=2))
            return sum(v for k, v in chain.items()
                       if all(f in a for f, a in zip(k, cochains))) % 2

        q, p = 3, 1
        cells = {d: list(combinations(range(6), d + 1)) for d in range(6)}
        # coboundary of the dual of {03, 04, 12, 14, 15, 34}
        alpha = {(0, 1, 2), (0, 1, 3), (0, 1, 5), (0, 2, 3), (0, 2, 4),
                 (0, 3, 4), (0, 3, 5), (0, 4, 5), (1, 2, 3), (1, 3, 5),
                 (2, 3, 4), (3, 4, 5)}
        top = (0, 1, 3, 4, 5)

        # left-hand side: Sq^1 Sq^1 at cochain level
        sq1 = {c for c in cells[3] if value(theta(1), (alpha, alpha), c)}
        sq1sq1 = {c for c in cells[4] if value(theta(2), (sq1, sq1), c)}
        lhs = int(top in sq1sq1)
        self.assertEqual(lhs, 1)

        coboundary = Surjection.adem_coboundary(q, p)
        rhs = 0
        for t in range(len(top)):
            face = top[:t] + top[t + 1:]
            rhs += value(coboundary, (alpha,) * 4, face)
        self.assertEqual(lhs, rhs % 2)

    def test_inputs(self):
        with self.assertRaises(ValueError):
            Steenrod.chain_operations(2, -1, -2, shape='not-a-shape')

        with self.assertRaises(ValueError):
            Steenrod.cartan_coboundaries(0, 3, shape='not-a-shape')

        with self.assertRaises(ValueError):
            Steenrod.adem_coboundaries(1, 0, 3, shape='not-a-shape')


if __name__ == '__main__':
    unittest.main()
