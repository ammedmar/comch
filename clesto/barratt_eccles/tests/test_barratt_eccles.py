import unittest
from clesto.barratt_eccles import BarrattEccles_element, BarrattEccles
from clesto.basics import SymmetricRing


class TestBarrattEccles_element(unittest.TestCase):
    def setUp(self):
        self.x = BarrattEccles_element({((1, 3, 2), (2, 3, 1)): 1})

    def test_arity(self):
        self.assertEqual(self.x.arity, 3)

    def test_degree(self):
        self.assertEqual(self.x.degree, 1)

    def test_boundary(self):
        self.assertEqual(self.x.boundary().boundary(), self.x.zero())

    def test_rmul(self):
        rho = SymmetricRing.rotation_element(3)
        a, b = (rho * self.x).boundary(), rho * self.x.boundary()
        self.assertEqual(a, b)

    def test_orbit(self):
        y = self.x.boundary().orbit(representation='trivial')
        z = self.x.orbit(representation='trivial').boundary()
        self.assertEqual(y, z.orbit('trivial'))

        y = self.x.boundary().orbit(representation='sign')
        z = self.x.orbit(representation='sign').boundary()
        self.assertEqual(y, z.orbit('sign'))

    def test_compose(self):
        i = 2
        x = BarrattEccles_element({((1, 2, 3), (3, 2, 1), (1, 3, 2)): 1})
        y = BarrattEccles_element({((1, 2, 3, 4), (4, 3, 2, 1), (4, 1, 3, 2)): 1})
        d_xy = x.compose(y, i).boundary()
        dx_y = x.boundary().compose(y, i)
        x_dy = x.compose(y.boundary(), i)
        self.assertEqual(d_xy, dx_y + (-1)**(x.degree) * x_dy)

    def test_table_reduction(self):
        b = BarrattEccles_element({((1, 2, 3, 4), (1, 4, 3, 2)): 1,
                                   ((1, 2, 4, 3), (3, 4, 2, 1)): 2})
        dtr_b = b.table_reduction().boundary()
        trd_b = b.boundary().table_reduction()
        self.assertEqual(dtr_b, trd_b)

    def test_alexander_whitney(self):
        pass


class TestBarrattEccles(unittest.TestCase):
    def test_steenrod_product(self):
        t = SymmetricRing.transposition_element(6)
        x = BarrattEccles.steenrod_product(6, 3).boundary()
        y = t * BarrattEccles.steenrod_product(6, 2)
        self.assertEqual(x, y)

        n = SymmetricRing.norm_element(5, torsion=7)
        x = BarrattEccles.steenrod_product(5, 6, torsion=7).boundary()
        y = n * BarrattEccles.steenrod_product(5, 5, torsion=7)
        self.assertEqual(x, y)


if __name__ == '__main__':
    unittest.main()
