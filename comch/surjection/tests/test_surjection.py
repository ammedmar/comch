import unittest
from math import comb

from comch.barratt_eccles import BarrattEccles
from comch.surjection import Surjection, SurjectionElement
from comch.symmetric import SymmetricRing, SymmetricRingElement
from comch.simplicial import SimplicialElement, Simplicial
from comch.cubical import Cubical


class TestSurjectionElement(unittest.TestCase):

    def setUp(self):
        self.x = SurjectionElement({(1, 3, 1, 2, 1): 1})

    def test_arity(self):
        self.assertEqual(self.x.arity, 3)

    def test_degree(self):
        self.assertEqual(self.x.degree, 2)

    def test_complexity(self):
        self.assertEqual(self.x.complexity, 2)

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

        self.assertEqual(3 * self.x, SurjectionElement({(1, 3, 1, 2, 1): 3}))

    def test_orbit(self):
        for conv in ['McClure-Smith', 'Berger-Fresse']:
            x = SurjectionElement({(1, 3, 1, 2, 1): 1}, convention=conv)
            a = x.orbit(representation='trivial')
            self.assertEqual(a, SurjectionElement({(1, 2, 1, 3, 1): 1},
                                                  convention=conv))

            b = x.orbit(representation='sign')
            self.assertEqual(b, SurjectionElement({(1, 2, 1, 3, 1): -1},
                                                  convention=conv))

    def test_call_simplicial(self):
        s = self.x
        s.convention = 'McClure-Smith'
        x = Simplicial.standard_element(3)
        ds_x = s.boundary()(x)
        d_sx = s(x).boundary()
        s_dx = s(x.boundary())
        self.assertEqual(d_sx - ((-1) ** s.degree) * s_dx, ds_x)

        x = SimplicialElement({((0, 1, 2), (3, 4), (5, 6)): 1})
        ds_x = s.boundary()(x, 2)
        d_sx = s(x, 2).boundary()
        s_dx = s(x.boundary(), 2)
        self.assertEqual(d_sx - ((-1) ** s.degree) * s_dx, ds_x)

    def test_boundary_berger_fresse_example(self):
        # worked example in Section 1.2.3 of Berger-Fresse
        u = SurjectionElement({(1, 3, 2, 1, 4, 2, 1): 1},
                              convention='Berger-Fresse')
        expected = SurjectionElement({
            (3, 2, 1, 4, 2, 1): 1,
            (1, 3, 1, 4, 2, 1): -1,
            (1, 3, 2, 4, 2, 1): 1,
            (1, 3, 2, 1, 4, 1): 1,
            (1, 3, 2, 1, 4, 2): -1,
        }, convention='Berger-Fresse')
        self.assertEqual(u.boundary(), expected)

    def test_compose_berger_fresse_example(self):
        # worked example in Section 1.2.6 of Berger-Fresse
        a = SurjectionElement({(1, 2, 1, 3): 1},
                              convention='Berger-Fresse')
        b = SurjectionElement({(1, 2, 1): 1},
                              convention='Berger-Fresse')
        expected = SurjectionElement({
            (1, 3, 1, 2, 1, 4): 1,
            (1, 2, 3, 2, 1, 4): -1,
            (1, 2, 1, 3, 1, 4): -1,
        }, convention='Berger-Fresse')
        self.assertEqual(a.compose(b, 1), expected)

    def test_compose_berger_fresse_derivation(self):
        x = SurjectionElement({(1, 2, 1, 3, 1): 1},
                              convention='Berger-Fresse')
        y = SurjectionElement({(2, 1, 2): 1},
                              convention='Berger-Fresse')
        for k in (1, 2, 3):
            d_xy = x.compose(y, k).boundary()
            dx_y = x.boundary().compose(y, k)
            x_dy = x.compose(y.boundary(), k)
            self.assertEqual(d_xy, dx_y + (-1) ** x.degree * x_dy)

    def test_call_simplicial_equivariance(self):
        # (sigma u)(x) = sigma u(x), permuting the output tensor factors
        # with the Koszul sign
        def permute_factors(element, sigma):
            answer = element.zero()
            for key, v in element.items():
                degrees = [len(spx) - 1 for spx in key]
                new_key = [None] * len(key)
                for idx, spx in enumerate(key):
                    new_key[sigma[idx] - 1] = spx
                sign_exp = sum(degrees[i] * degrees[j]
                               for i in range(len(key))
                               for j in range(i + 1, len(key))
                               if sigma[i] > sigma[j])
                answer += answer.create(
                    {tuple(new_key): (-1) ** (sign_exp % 2) * v})
            return answer

        x = Simplicial.standard_element(3)
        for convention in ('Berger-Fresse', 'McClure-Smith'):
            u = SurjectionElement({(1, 2, 3, 1, 2): 1},
                                  convention=convention)
            for sigma in ((2, 3, 1), (3, 1, 2), (2, 1, 3)):
                su = SymmetricRingElement({sigma: 1}) * u
                self.assertEqual(su(x), permute_factors(u(x), sigma))

    def test_call_simplicial_naturality(self):
        # the action commutes with the simplicial map relabeling vertices
        u = SurjectionElement({(1, 2, 1, 3): 1},
                              convention='Berger-Fresse')
        labels = (0, 2, 5, 6)
        pushed = SimplicialElement()
        for k, v in u(Simplicial.standard_element(3)).items():
            new_k = tuple(tuple(labels[i] for i in spx) for spx in k)
            pushed += pushed.create({new_k: v})
        self.assertEqual(pushed, u(SimplicialElement({(labels,): 1})))

    def test_call_simplicial_iterated_diagonal(self):
        # AW(1, 2, ..., r) is the (r-1)-fold Alexander-Whitney diagonal
        x = Simplicial.standard_element(3)
        for r in (2, 3, 4):
            u = SurjectionElement({tuple(range(1, r + 1)): 1},
                                  convention='Berger-Fresse')
            self.assertEqual(u(x), x.iterated_diagonal(r - 1))

    def test_call_simplicial_berger_fresse(self):
        s = SurjectionElement({(1, 3, 1, 2, 1): 1},
                              convention='Berger-Fresse')
        x = Simplicial.standard_element(3)
        ds_x = s.boundary()(x)
        d_sx = s(x).boundary()
        s_dx = s(x.boundary())
        self.assertEqual(d_sx - ((-1) ** s.degree) * s_dx, ds_x)

        x = SimplicialElement({((0, 1, 2), (3, 4), (5, 6)): 1})
        ds_x = s.boundary()(x, 2)
        d_sx = s(x, 2).boundary()
        s_dx = s(x.boundary(), 2)
        self.assertEqual(d_sx - ((-1) ** s.degree) * s_dx, ds_x)

    def test_call_simplicial_berger_fresse_cup_products(self):
        # closed formulas at the end of Section 2.2 of Berger-Fresse
        n = 4
        x = Simplicial.standard_element(n)

        cup_1 = SurjectionElement({(1, 2, 1): 1},
                                  convention='Berger-Fresse')
        expected = SimplicialElement()
        for i in range(n + 1):
            for j in range(i + 1, n + 1):
                key = (tuple(range(i + 1)) + tuple(range(j, n + 1)),
                       tuple(range(i, j + 1)))
                sign = (-1) ** (((n - j) * (j - i) + i) % 2)
                expected += expected.create({key: sign})
        self.assertEqual(cup_1(x), expected)

        cup_2 = SurjectionElement({(1, 2, 1, 2): 1},
                                  convention='Berger-Fresse')
        expected = SimplicialElement()
        for i in range(n + 1):
            for j in range(i + 1, n + 1):
                for k in range(j + 1, n + 1):
                    key = (tuple(range(i + 1)) + tuple(range(j, k + 1)),
                           tuple(range(i, j + 1)) + tuple(range(k, n + 1)))
                    sign = (-1) ** (((k - j) * (j - i + 1) + i + j) % 2)
                    expected += expected.create({key: sign})
        self.assertEqual(cup_2(x), expected)

    def test_call_simplicial_conventions_agree_mod_two(self):
        for key in [(1, 2, 1), (2, 1, 3, 1), (1, 2, 3, 1, 2)]:
            bf = SurjectionElement({key: 1}, torsion=2,
                                   convention='Berger-Fresse')
            ms = SurjectionElement({key: 1}, torsion=2,
                                   convention='McClure-Smith')
            x = Simplicial.standard_element(3, torsion=2)
            self.assertEqual(bf(x), ms(x))

    def test_call_simplicial_berger_fresse_operadic(self):
        # (u o_k v)(x) = (-1)^{|u||v|} v(u(x), k), the Koszul sign of the
        # coendomorphism operad composition
        u = SurjectionElement({(2, 1, 2): 1}, convention='Berger-Fresse')
        v = SurjectionElement({(1, 2, 1): 1}, convention='Berger-Fresse')
        x = Simplicial.standard_element(3)
        for k in (1, 2):
            sign = (-1) ** (u.degree * v.degree)
            self.assertEqual(u.compose(v, k)(x), sign * v(u(x), k))

    def test_call_cubical(self):
        s = self.x
        y = Cubical.standard_element(3)
        ds_y = s.boundary()(y)
        d_sy = s(y).boundary()
        sdy = s(y.boundary())
        self.assertEqual(d_sy - ((-1) ** s.degree) * sdy, ds_y)

        y = Cubical.standard_element(3, 2)
        ds_y = s.boundary()(y, 2)
        d_sy = s(y, 2).boundary()
        sdy = s(y.boundary(), 2)
        self.assertEqual(d_sy - ((-1) ** s.degree) * sdy, ds_y)

    def test_compose_bf(self):
        i = 3
        x = SurjectionElement({(3, 2, 1, 2, 1, 3): 1},
                              convention='Berger-Fresse')
        y = SurjectionElement({(3, 1, 2, 1, 4, 3): 1},
                              convention='Berger-Fresse')
        dx, dy = x.boundary(), y.boundary()
        dx_y, x_dy = dx.compose(y, i), x.compose(dy, i)
        xy = x.compose(y, i)
        d_xy = xy.boundary()
        self.assertEqual(d_xy - dx_y - (-1) ** (x.degree) * x_dy, x.zero())

    def test_suspension(self):
        x = SurjectionElement({(1, 3, 2, 1, 2, 3, 4): 1},
                              convention='Berger-Fresse')
        y = SurjectionElement({(1, 2, 3, 1, 2, 3): 1},
                              convention='Berger-Fresse')
        sx = x.suspension()
        sy = y.suspension()
        xy = x.compose(y, 1)
        self.assertEqual(xy.suspension(), sx.compose(sy, 1))

        z = SurjectionElement({(1, 2, 1, 3, 1, 2, 3): 1},
                              convention='Berger-Fresse')
        sz = z.suspension()
        xz = x.compose(z, 2)
        self.assertEqual(xz.suspension(), sx.compose(sz, 1))


class TestSurjection(unittest.TestCase):
    def test_deprecated_aliases(self):
        with self.assertWarns(DeprecationWarning):
            old_map = Surjection.may_steenrod_structure(
                3, 2, torsion=3, convention='Berger-Fresse')
        self.assertEqual(
            old_map,
            Surjection.may_steenrod_map(
                3, 2, torsion=3, convention='Berger-Fresse'),
        )

        with self.assertWarns(DeprecationWarning):
            old_operation = Surjection.steenrod_operation(
                3, 2, 3, convention='Berger-Fresse')
        self.assertEqual(
            old_operation,
            Surjection.steenrod_operation_representative(
                3, 2, 3, convention='Berger-Fresse'),
        )

    def test_may_steenrod_map(self):
        arity = 6
        t = SymmetricRing.transposition_element(arity)
        x = Surjection.may_steenrod_map(arity, 3).boundary()
        y = t * Surjection.may_steenrod_map(arity, 2)
        self.assertEqual(x, y)

        arity = 5
        n = SymmetricRing.norm_element(arity, torsion=7)
        x = Surjection.may_steenrod_map(arity, 4, torsion=7).boundary()
        y = n * Surjection.may_steenrod_map(arity, 3, torsion=7)
        self.assertEqual(x, y)

        t = SymmetricRing.transposition_element(arity)
        x = Surjection.may_steenrod_map(
            arity, 3, convention='McClure-Smith').boundary()
        y = t * Surjection.may_steenrod_map(
            arity, 2, convention='McClure-Smith')
        self.assertEqual(x, y)

        self.assertEqual(Surjection.may_steenrod_map(2, 1, torsion=0).torsion,
                         0)
        self.assertFalse(Surjection.may_steenrod_map(2, 1, torsion=1))
        self.assertEqual(Surjection.may_steenrod_map(2, 1, torsion=-5),
                         Surjection.may_steenrod_map(2, 1, torsion=5))
        with self.assertRaises(TypeError):
            Surjection.may_steenrod_map(2, 1, torsion=False)

    def test_steenrod_operation_representative(self):
        x = Surjection.steenrod_operation_representative(3, 2, 3)
        y = SurjectionElement({
            (1, 2, 3, 1, 2): 2,
            (1, 3, 1, 2, 3): 2,
            (1, 2, 3, 2, 3): 2,
        }, torsion=3, convention='McClure-Smith')
        self.assertEqual(x, y)

        x = Surjection.steenrod_operation_representative(
            2, 3, 2, convention='Berger-Fresse')
        self.assertEqual(x.convention, 'Berger-Fresse')

        self.assertEqual(
            Surjection.steenrod_operation_representative(2, 0, 1),
            SurjectionElement(torsion=2, convention='McClure-Smith'),
        )
        with self.assertRaises(TypeError):
            Surjection.steenrod_operation_representative(2, 0, 0,
                                                         bockstein=True)
        with self.assertRaises(TypeError):
            Surjection.steenrod_operation_representative(2, 0, 0,
                                                         bockstein=1)

    def test_cartan_homotopy(self):
        expected_0 = SurjectionElement({
            (1, 3, 2, 3, 4): 1,
        }, torsion=2, convention='Berger-Fresse')
        self.assertEqual(Surjection.cartan_homotopy(0), expected_0)

        for convention in ('Berger-Fresse', 'McClure-Smith'):
            for i in range(4):
                answer = Surjection.cartan_homotopy(i, convention=convention)
                expected = BarrattEccles.cartan_homotopy(i).table_reduction()
                expected.convention = convention
                self.assertEqual(answer, expected)
                self.assertEqual(answer.convention, convention)

        with self.assertRaises(TypeError):
            Surjection.cartan_homotopy(1.0)
        with self.assertRaises(TypeError):
            Surjection.cartan_homotopy(True)
        with self.assertRaises(ValueError):
            Surjection.cartan_homotopy(-1)
        with self.assertRaises(ValueError):
            Surjection.cartan_homotopy(0, convention='unknown')

        odd = Surjection.cartan_homotopy(
            2, convention='Berger-Fresse', arity=3)
        self.assertEqual(
            odd,
            BarrattEccles.cartan_homotopy(
                2, arity=3).table_reduction(),
        )
        self.assertEqual(odd.torsion, 3)
        self.assertEqual(odd.arity, 6)
        self.assertEqual(odd.degree, 3)
        with self.assertRaises(NotImplementedError):
            Surjection.cartan_homotopy(
                2, convention='McClure-Smith', arity=3)

    def test_adem_homotopy(self):
        expected = SurjectionElement({
            (1, 2, 3, 2, 4, 2): 1,
            (1, 3, 4, 1, 4, 2): 1,
        }, torsion=2, convention='Berger-Fresse')
        self.assertEqual(Surjection.adem_homotopy(1, 0), expected)

        for convention in ('Berger-Fresse', 'McClure-Smith'):
            for q, p in ((1, 0), (2, 1)):
                answer = Surjection.adem_homotopy(
                    q, p, convention=convention)
                expected = BarrattEccles.adem_homotopy(q, p).table_reduction()
                expected.convention = convention
                self.assertEqual(answer, expected)
                self.assertEqual(answer.convention, convention)

        with self.assertRaises(TypeError):
            Surjection.adem_homotopy(1.0, 0)
        with self.assertRaises(TypeError):
            Surjection.adem_homotopy(True, 0)
        with self.assertRaises(ValueError):
            Surjection.adem_homotopy(0, 0)
        with self.assertRaises(ValueError):
            Surjection.adem_homotopy(1, -1)
        with self.assertRaises(ValueError):
            Surjection.adem_homotopy(1, 0, convention='unknown')

    def test_adem_coboundary(self):
        expected = SurjectionElement({
            (1, 2, 3, 2, 4, 2): 1,
            (1, 3, 4, 1, 4, 2): 1,
            (1, 3, 4, 3, 1, 2): 1,
            (1, 2, 3, 4, 3, 2): 1,
        }, torsion=2, convention='Berger-Fresse')
        self.assertEqual(Surjection.adem_coboundary(1, 0), expected)

        for convention in ('Berger-Fresse', 'McClure-Smith'):
            for q, p in ((1, 0), (2, 1)):
                answer = Surjection.adem_coboundary(
                    q, p, convention=convention)
                expected = BarrattEccles.adem_coboundary(
                    q, p).table_reduction()
                expected.convention = convention
                self.assertEqual(answer, expected)
                self.assertEqual(answer.convention, convention)

        with self.assertRaises(TypeError):
            Surjection.adem_coboundary(1.0, 0)
        with self.assertRaises(TypeError):
            Surjection.adem_coboundary(True, 0)
        with self.assertRaises(ValueError):
            Surjection.adem_coboundary(0, 0)
        with self.assertRaises(ValueError):
            Surjection.adem_coboundary(1, -1)
        with self.assertRaises(ValueError):
            Surjection.adem_coboundary(1, 0, convention='unknown')

    def test_adem_coboundary_identity(self):
        r"""The complete element has the boundary asserted in [BMM].

        The equality is in coinvariants because evaluation on
        alpha^{otimes 4} is invariant under the Sigma_4 action.
        """

        def theta(degree):
            return Surjection.may_steenrod_map(
                2, degree, torsion=2, convention='Berger-Fresse')

        def symmetric_adem_term(outer_index, inner_index):
            answer = SurjectionElement(
                torsion=2, convention='Berger-Fresse')
            for ell in range(inner_index // 2 + 1):
                coefficient = comb(
                    inner_index - ell, inner_index - 2 * ell) % 2
                outer_degree = outer_index - inner_index + 2 * ell
                if coefficient == 0 or outer_degree < 0:
                    continue
                inner = theta(inner_index - ell)
                answer += theta(outer_degree).compose(
                    inner, 1).compose(inner, 3)
            return answer

        for q, p in ((1, 0), (2, 1), (3, 1)):
            coboundary = Surjection.adem_coboundary(
                q, p, convention='Berger-Fresse')
            residual = coboundary.boundary()
            residual += symmetric_adem_term(q, p)
            residual += symmetric_adem_term(p, q)
            self.assertFalse(residual.orbit(), msg=(q, p))


if __name__ == '__main__':
    unittest.main()
