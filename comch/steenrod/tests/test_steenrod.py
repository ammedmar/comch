import unittest
from itertools import combinations

from comch.cubical import CubicalElement
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
        }, torsion=2)
        self.assertEqual(x, y)

        self.assertEqual(Steenrod.adem_coboundaries(1, 0, 3),
                         SimplicialElement(torsion=2))

        z = Steenrod.adem_coboundaries(1, 0, 1, shape='cube')
        self.assertIsInstance(z, CubicalElement)
        self.assertEqual(len(z), 0)

    def test_inputs(self):
        with self.assertRaises(ValueError):
            Steenrod.chain_operations(2, -1, -2, shape='not-a-shape')

        with self.assertRaises(ValueError):
            Steenrod.cartan_coboundaries(0, 3, shape='not-a-shape')

        with self.assertRaises(ValueError):
            Steenrod.adem_coboundaries(1, 0, 3, shape='not-a-shape')


if __name__ == '__main__':
    unittest.main()
