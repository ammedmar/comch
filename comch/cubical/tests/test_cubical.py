import unittest
from comch.cubical import CubicalElement, Cubical
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
        rho = SymmetricRing.rotation_element(2)
        a, b = (rho * self.x).boundary(), rho * self.x.boundary()
        self.assertEqual(a, b)

    def test_iterated_diagonal(self):
        delta_d_x = self.x.boundary().iterated_diagonal(3, 2)
        d_delta_x = self.x.iterated_diagonal(3, 2).boundary()
        self.assertEqual(delta_d_x, d_delta_x)

    def test_join(self):
        x, y = (2, 0), (0, 1)
        a = CubicalElement({(x,): 1})
        da = a.boundary()
        b = CubicalElement({(y,): 1})
        db = b.boundary()
        c = CubicalElement({(x, y): 1})
        da_b = CubicalElement({(p[0], y): r for p, r in da.items()})
        a_db = CubicalElement({(x, q[0]): r for q, r in db.items()})
        left = c.join().boundary() + da_b.join() + (-1) ** b.degree * a_db.join()
        self.assertEqual(left, - a)


if __name__ == '__main__':
    unittest.main()
