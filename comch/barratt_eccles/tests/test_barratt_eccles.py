import unittest
from comch.barratt_eccles import BarrattEccles, BarrattEcclesElement
from comch.symmetric import SymmetricRing, SymmetricRingElement


class TestBarrattEcclesElement(unittest.TestCase):
    def setUp(self):
        self.x = BarrattEcclesElement({((1, 3, 2), (2, 3, 1)): 1})

    def test_arity(self):
        self.assertEqual(self.x.arity, 3)

    def test_complexity(self):
        self.assertEqual(self.x.complexity, 1)

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
        x = BarrattEcclesElement({((1, 2, 3), (3, 2, 1), (1, 3, 2)): 1})
        y = BarrattEcclesElement({((1, 2, 3, 4), (4, 3, 2, 1), (4, 1, 3, 2)): 1})
        d_xy = x.compose(y, i).boundary()
        dx_y = x.boundary().compose(y, i)
        x_dy = x.compose(y.boundary(), i)
        self.assertEqual(d_xy, dx_y + (-1)**(x.degree) * x_dy)

    def test_table_reduction(self):
        # chain map
        b = BarrattEcclesElement({((1, 2, 3, 4), (1, 4, 3, 2)): 1,
                                  ((1, 2, 4, 3), (3, 4, 2, 1)): 2})
        dtr_b = b.table_reduction().boundary()
        trd_b = b.boundary().table_reduction()
        self.assertEqual(dtr_b, trd_b)

        # equivariant
        sym = SymmetricRingElement({(1, 3, 2, 4): 1,
                                    (2, 1, 3, 4): -1})
        sym_tr_b = sym * (b.table_reduction())
        tr_sym_b = (sym * b).table_reduction()
        self.assertEqual(sym_tr_b, tr_sym_b)

        a = BarrattEcclesElement({((1, 2, 3), (1, 3, 2)): 1,
                                  ((2, 1, 3), (1, 2, 3)): 2})
        tr_ab = a.compose(b, 2).table_reduction()
        tra_trb = a.table_reduction().compose(b.table_reduction(), 2)
        self.assertEqual(tr_ab, tra_trb)

    def test_alexander_whitney(self):
        pass


class TestBarrattEccles(unittest.TestCase):
    def test_cartan_homotopy(self):
        identity = (1, 2, 3, 4)
        transposition_23 = (1, 3, 2, 4)
        transposition_34 = (1, 2, 4, 3)
        double_transposition = (2, 1, 4, 3)
        conjugated_transposition = (2, 4, 1, 3)

        expected_0 = BarrattEcclesElement({
            (transposition_23, identity): 1,
        }, torsion=2)
        self.assertEqual(BarrattEccles.cartan_homotopy(0), expected_0)

        expected_1 = BarrattEcclesElement({
            (transposition_23, identity, double_transposition): 1,
            (transposition_23, conjugated_transposition,
             double_transposition): 1,
            (identity, transposition_34, double_transposition): 1,
        }, torsion=2)
        self.assertEqual(BarrattEccles.cartan_homotopy(1), expected_1)

        expected_table_reduction_2 = {
            (1, 3, 2, 3, 4, 3, 4),
            (1, 2, 4, 1, 4, 3, 4),
            (1, 2, 1, 3, 2, 3, 4),
            (1, 3, 2, 3, 2, 3, 4),
        }
        self.assertEqual(
            set(BarrattEccles.cartan_homotopy(2).table_reduction()),
            expected_table_reduction_2,
        )

        with self.assertRaises(TypeError):
            BarrattEccles.cartan_homotopy(1.0)
        with self.assertRaises(TypeError):
            BarrattEccles.cartan_homotopy(True)
        with self.assertRaises(ValueError):
            BarrattEccles.cartan_homotopy(-1)

    def test_cartan_homotopy_boundary(self):
        x_0 = BarrattEcclesElement({((1, 2),): 1}, torsion=2)
        transposition_23 = SymmetricRingElement(
            {(1, 3, 2, 4): 1}, torsion=2)
        double_transposition_plus_identity = SymmetricRingElement({
            (2, 1, 4, 3): 1,
            (1, 2, 3, 4): 1,
        }, torsion=2)

        for i in range(5):
            x_i = BarrattEccles.may_steenrod_map(2, i, torsion=2)
            first = transposition_23 * x_i.compose(x_0, 1).compose(x_0, 3)

            second = BarrattEcclesElement(torsion=2)
            for (left, right), coefficient in x_i.diagonal().items():
                left = BarrattEcclesElement({left: 1}, torsion=2)
                right = BarrattEcclesElement({right: 1}, torsion=2)
                second += coefficient * x_0.compose(left, 1).compose(right, 3)

            boundary = BarrattEccles.cartan_homotopy(i).boundary()
            if i > 0:
                boundary += double_transposition_plus_identity * (
                    BarrattEccles.cartan_homotopy(i - 1))
            self.assertEqual(boundary, first + second)

    def test_may_steenrod_map(self):
        t = SymmetricRing.transposition_element(6)
        x = BarrattEccles.may_steenrod_map(6, 3).boundary()
        y = t * BarrattEccles.may_steenrod_map(6, 2)
        self.assertEqual(x, y)

        n = SymmetricRing.norm_element(5, torsion=7)
        x = BarrattEccles.may_steenrod_map(5, 6, torsion=7).boundary()
        y = n * BarrattEccles.may_steenrod_map(5, 5, torsion=7)
        self.assertEqual(x, y)

        self.assertEqual(BarrattEccles.may_steenrod_map(2, 1, torsion=0).torsion,
                         0)
        self.assertFalse(BarrattEccles.may_steenrod_map(2, 1, torsion=1))
        self.assertEqual(BarrattEccles.may_steenrod_map(2, 1, torsion=-5),
                         BarrattEccles.may_steenrod_map(2, 1, torsion=5))
        with self.assertRaises(TypeError):
            BarrattEccles.may_steenrod_map(2, 1, torsion=False)

    def test_steenrod_operation_representative(self):
        x = BarrattEccles.steenrod_operation_representative(3, 2, 3)
        y = BarrattEcclesElement({
            ((1, 2, 3), (2, 3, 1), (3, 1, 2)): 2,
            ((1, 2, 3), (3, 1, 2), (1, 2, 3)): 2,
        }, torsion=3)
        self.assertEqual(x, y)

        self.assertEqual(
            BarrattEccles.steenrod_operation_representative(2, 0, 1),
            BarrattEcclesElement(torsion=2),
        )
        with self.assertRaises(TypeError):
            BarrattEccles.steenrod_operation_representative(
                2, 0, 0, bockstein=True)
        with self.assertRaises(TypeError):
            BarrattEccles.steenrod_operation_representative(
                2, 0, 0, bockstein=1)
        with self.assertRaises(TypeError):
            BarrattEccles.steenrod_operation_representative(
                2, 0, 0, convention='McClure-Smith')


if __name__ == '__main__':
    unittest.main()
