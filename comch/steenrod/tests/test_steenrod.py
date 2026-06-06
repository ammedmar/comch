import unittest

from comch.cubical import CubicalElement
from comch.simplicial import SimplicialElement
from comch.steenrod import Steenrod


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

        z = Steenrod.cartan_coboundaries(0, 3, shape='cube')
        self.assertIsInstance(z, CubicalElement)
        self.assertEqual(len(z), 6)

    def test_inputs(self):
        with self.assertRaises(ValueError):
            Steenrod.chain_operations(2, -1, -2, shape='not-a-shape')

        with self.assertRaises(ValueError):
            Steenrod.cartan_coboundaries(0, 3, shape='not-a-shape')


if __name__ == '__main__':
    unittest.main()
