import unittest
from clesto.surjection import Surjection_element
from clesto.basics import SymmetricRing
from clesto.eilenberg_zilber import EilenbergZilber, CubicalEilenbergZilber


class TestSurjection_element(unittest.TestCase):

    def setUp(self):
        self.x = Surjection_element({(1, 3, 1, 2, 1): 1})

    def test_arity(self):
        self.assertEqual(self.x.arity, 3)

    def test_degree(self):
        self.assertEqual(self.x.degree, 2)

    def test_complexity(self):
        self.assertEqual(self.x.complexity, 1)

    def test_boundary(self):
        self.assertEqual(self.x.boundary().boundary(), self.x.zero())

        self.x.convention = 'Berger-Fresse'
        self.assertEqual(self.x.boundary().boundary(), self.x.zero())

        self.x.set_torsion(2)
        self.assertEqual(self.x.boundary().boundary(), self.x.zero())

    def test_rmul(self):
        rho = SymmetricRing.rotation_element(3)
        a, b = (rho * self.x).boundary(), rho * self.x.boundary()
        self.assertEqual(a, b)

        self.x.convention = 'Berger-Fresse'
        a, b = (rho * self.x).boundary(), rho * self.x.boundary()
        self.assertEqual(a, b)

    def test_orbit(self):
        a = self.x.orbit(representation='trivial')
        self.assertEqual(a, Surjection_element({(1, 2, 1, 3, 1): 1}))

        b = self.x.orbit(representation='sign')
        self.assertEqual(b, Surjection_element({(1, 2, 1, 3, 1): -1}))

    def test_call(self):
        s = self.x
        x = EilenbergZilber.standard_element(3)
        ds_x = s.boundary()(x)
        d_sx = s(x).boundary()
        sdx = s(x.boundary())
        self.assertEqual(d_sx - ((-1)**(s.degree)) * sdx, ds_x)

        y = CubicalEilenbergZilber.standard_element(3)
        ds_y = s.boundary()(y)
        d_sy = s(y).boundary()
        sdy = s(y.boundary())
        self.assertEqual(d_sy - ((-1)**(s.degree)) * sdy, ds_y)

    def test_compose_bf(self):
        i = 3
        x = Surjection_element({(3, 2, 1, 2, 1, 3): 1}, convention='Berger-Fresse')
        y = Surjection_element({(3, 1, 2, 1, 4, 3): 1}, convention='Berger-Fresse')
        dx = x.boundary()
        dy = y.boundary()
        xy = x.compose(y, i)
        d_xy = xy.boundary()
        dx_y = dx.compose(y, i)
        x_dy = x.compose(dy, i)
        self.assertEqual(d_xy - dx_y - (-1)**(x.degree) * x_dy, x.zero())


if __name__ == '__main__':
    unittest.main()
