import unittest

from comch.simplicial import SimplicialElement
from comch.steenrod import Steenrod


class TestSteenrod(unittest.TestCase):
    def test_chain(self):
        x = Steenrod.chain(2, -1, -2)
        y = SimplicialElement({
            ((0, 2, 3), (0, 1, 2)): 1,
            ((0, 1, 3), (1, 2, 3)): 1,
        }, torsion=2)
        self.assertEqual(x, y)

    def test_inputs(self):
        with self.assertRaises(ValueError):
            Steenrod.chain(2, -1, -2, shape='not-a-shape')


if __name__ == '__main__':
    unittest.main()
