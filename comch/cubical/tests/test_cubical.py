import unittest
from comch.cubical import CubicalElement
from comch.symmetric import SymmetricRing


class TestBarrattEcclesElement(unittest.TestCase):
    def setUp(self):
        self.x = CubicalElement({((0, 2, 1), (0, 1, 2)): 1,
                                 ((2, 1, 1), (1, 2, 0)): -1,
                                 ((0, 1, 0), (2, 2, 1)): 1})

    def test_arity(self):
        self.assertEqual(self.x.arity, 2)

    def test_degree(self):
        self.assertEqual(self.x.degree, 2)

    def test_boundary(self):
        self.assertEqual(self.x.boundary().boundary(), self.x.zero())

    def test_rmul(self):
        pass

    def test_iterated_diagonal(self):
        delta_d_x = self.x.boundary().iterated_diagonal(3, 2)
        d_delta_x = self.x.iterated_diagonal(3, 2).boundary()
        self.assertEqual(delta_d_x, d_delta_x)


if __name__ == '__main__':
    unittest.main()
