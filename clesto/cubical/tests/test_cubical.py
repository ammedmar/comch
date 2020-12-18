import unittest
from clesto.simplicial import SimplicialEZ_element
from clesto.symmetric import SymmetricRing


class TestBarrattEccles_element(unittest.TestCase):
    def setUp(self):
        self.x = SimplicialEZ_element({((0,), (0, 1, 2)): 1,
                                       ((0, 1), (1, 2)): -1,
                                       ((0, 1, 2), (2,)): 1})

    def test_arity(self):
        self.assertEqual(self.x.arity, 2)

    def test_degree(self):
        self.assertEqual(self.x.degree, 2)

    def test_boundary(self):
        self.assertEqual(self.x.boundary().boundary(), self.x.zero())

    def test_rmul(self):
        rho = SymmetricRing.rotation_element(2)
        a, b = (rho * self.x).boundary(), rho * self.x.boundary()
        self.assertEqual(a, b)

    def test_iterated_diagonal(self):
        pass


if __name__ == '__main__':
    unittest.main()
