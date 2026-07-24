import unittest
from math import comb

from comch.barratt_eccles import BarrattEccles, BarrattEcclesElement
from comch.cubical import Cubical
from comch.simplicial import Simplicial, SimplicialElement
from comch.symmetric import (
    SymmetricGroupElement, SymmetricRing, SymmetricRingElement)


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

    def test_table_reduction_filtration(self):
        # TR preserves the little-cubes filtration ([BF] Section 1.6): an
        # element of complexity n - 1 lies in F_n E and its table reduction
        # lies in F_n X, formed by surjections of complexity at most n. The
        # anchor examples are those of [BF]: ((1,2),(2,1),(1,2)) in F_3 E
        # and its table reduction (1,2,1,2) in F_3 X.
        b = BarrattEcclesElement({((1, 2), (2, 1), (1, 2)): 1})
        self.assertEqual(b.complexity, 2)
        self.assertEqual(b.table_reduction().complexity, 3)

        keys = [((1, 2, 3), (2, 1, 3), (3, 1, 2), (3, 2, 1)),
                ((2, 3, 1), (3, 1, 2)),
                ((1, 2, 3, 4), (1, 4, 3, 2), (1, 2, 4, 3))]
        for key in keys:
            b = BarrattEcclesElement({key: 1})
            tr = b.table_reduction()
            if tr:
                self.assertLessEqual(tr.complexity, b.complexity + 1)

    def test_call_simplicial_chain_map(self):
        x = BarrattEcclesElement({((1, 2, 3), (1, 3, 2), (2, 3, 1)): 1})
        c = Simplicial.standard_element(3)
        dx_c = x.boundary()(c)
        d_xc = x(c).boundary()
        x_dc = x(c.boundary())
        self.assertEqual(d_xc - ((-1) ** x.degree) * x_dc, dx_c)

    def test_call_simplicial_factors_through_table_reduction(self):
        x = BarrattEcclesElement({((1, 2), (2, 1), (1, 2)): 1})
        c = SimplicialElement({((0, 1, 2), (3, 4, 5)): 1})
        for coord in (1, 2):
            self.assertEqual(x(c, coord),
                             x.table_reduction()(c, coord))

    def test_call_cubical_not_implemented(self):
        x = BarrattEcclesElement({((1, 2), (2, 1)): 1})
        with self.assertRaises(NotImplementedError):
            x(Cubical.standard_element(2))


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
        with self.assertRaises(TypeError):
            BarrattEccles.cartan_homotopy(0, arity=True)
        with self.assertRaises(ValueError):
            BarrattEccles.cartan_homotopy(0, arity=1)

        odd = BarrattEccles.cartan_homotopy(2, arity=3)
        self.assertEqual(odd.torsion, 3)
        self.assertEqual(odd.arity, 6)
        self.assertEqual(odd.degree, 3)
        self.assertEqual(len(odd), 16)
        self.assertEqual(len(odd.table_reduction()), 30)

    def test_odd_cartan_relator_identity(self):
        arity = 3
        rotation = SymmetricGroupElement((2, 3, 1))
        block_rotation = SymmetricGroupElement((2, 3, 1, 5, 6, 4))
        tau = SymmetricGroupElement((1, 3, 5, 2, 4, 6))
        product = BarrattEccles.may_steenrod_map(
            2, 0, torsion=arity)

        def right_action(element, permutation, exponent=1):
            return element.create({
                tuple(value * (permutation ** exponent)
                      for value in simplex): coefficient
                for simplex, coefficient in element.items()
            })

        def interleaved_product(element):
            answer = element
            for position in range(1, 2 * arity, 2):
                answer = answer.compose(product, position)
            return answer

        def reordered_product(degree):
            source = BarrattEccles.may_steenrod_map(
                arity, degree, torsion=arity)
            answer = BarrattEcclesElement(torsion=arity)
            for simplex, coefficient in source.items():
                image = []
                for permutation in simplex:
                    vertex = interleaved_product(BarrattEcclesElement({
                        (permutation,): 1}, torsion=arity))
                    image.append(next(iter(vertex))[0] * tau)
                answer += answer.create({tuple(image): coefficient})
            return answer

        def cyclic_diagonal_product(degree):
            half_degree = degree // 2
            terms = []
            if degree % 2 == 0:
                terms.extend(
                    (2 * j, 2 * (half_degree - j), 0, 0)
                    for j in range(half_degree + 1)
                )
                terms.extend(
                    (2 * j + 1, 2 * (half_degree - 1 - j) + 1,
                     left_action, right_action_power)
                    for j in range(half_degree)
                    for left_action in range(arity)
                    for right_action_power in range(
                        left_action + 1, arity)
                )
            else:
                for j in range(half_degree + 1):
                    terms.extend((
                        (2 * j, 2 * (half_degree - j) + 1, 0, 0),
                        (2 * j + 1, 2 * (half_degree - j), 0, 1),
                    ))

            answer = BarrattEcclesElement(torsion=arity)
            for left_degree, right_degree, left_power, right_power in terms:
                left = right_action(
                    BarrattEccles.may_steenrod_map(
                        arity, left_degree, torsion=arity),
                    rotation, left_power)
                right = right_action(
                    BarrattEccles.may_steenrod_map(
                        arity, right_degree, torsion=arity),
                    rotation, right_power)
                answer += product.compose(left, 1).compose(
                    right, arity + 1)
            return answer

        for degree in range(4):
            homotopy = BarrattEccles.cartan_homotopy(
                degree, arity=arity)
            left_hand_side = homotopy.boundary()
            if degree:
                previous = BarrattEccles.cartan_homotopy(
                    degree - 1, arity=arity)
                if degree % 2:
                    left_hand_side += (
                        right_action(previous, block_rotation) - previous)
                else:
                    for exponent in range(arity):
                        left_hand_side += right_action(
                            previous, block_rotation, exponent)

            right_hand_side = (
                cyclic_diagonal_product(degree)
                - reordered_product(degree)
            )
            self.assertEqual(left_hand_side, right_hand_side)

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

    def test_adem_homotopy(self):
        identity = (1, 2, 3, 4)
        transposition_23 = (1, 3, 2, 4)
        transposition_34 = (1, 2, 4, 3)
        double_transposition = (2, 1, 4, 3)
        d8_rotation = (3, 4, 1, 2)
        transposition_13_24_34 = (3, 1, 4, 2)

        expected = BarrattEcclesElement({
            (identity, transposition_34, double_transposition): 1,
            (identity, transposition_23, transposition_13_24_34): 1,
            (identity, d8_rotation, transposition_13_24_34): 1,
        }, torsion=2)
        self.assertEqual(BarrattEccles.adem_homotopy(1, 0), expected)

        expected_table_reduction = {
            (1, 2, 3, 2, 4, 2),
            (1, 3, 4, 1, 4, 2),
        }
        self.assertEqual(
            set(BarrattEccles.adem_homotopy(1, 0).table_reduction()),
            expected_table_reduction,
        )

        for q, p in ((1, 0), (2, 1)):
            homotopy = BarrattEccles.adem_homotopy(q, p)
            self.assertEqual(homotopy.arity, 4)
            self.assertEqual(homotopy.degree, p + q + 1)
            self.assertEqual(homotopy.boundary().boundary(), homotopy.zero())
            self.assertEqual(
                homotopy.table_reduction().boundary(),
                homotopy.boundary().table_reduction(),
            )

        with self.assertRaises(TypeError):
            BarrattEccles.adem_homotopy(1.0, 0)
        with self.assertRaises(TypeError):
            BarrattEccles.adem_homotopy(True, 0)
        with self.assertRaises(ValueError):
            BarrattEccles.adem_homotopy(0, 0)
        with self.assertRaises(ValueError):
            BarrattEccles.adem_homotopy(1, -1)

    def test_adem_coboundary(self):
        def theta(degree):
            return BarrattEccles.may_steenrod_map(
                2, degree, torsion=2)

        def symmetric_adem_term(outer_index, inner_index):
            answer = BarrattEcclesElement(torsion=2)
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

        homotopy = BarrattEccles.adem_homotopy(1, 0)
        theta_0 = theta(0)
        theta_1 = theta(1)
        correction = theta_1.compose(theta_0, 1).compose(theta_1, 3)
        self.assertEqual(
            BarrattEccles.adem_coboundary(1, 0),
            homotopy + correction,
        )

        for q, p in ((1, 0), (2, 1), (3, 1)):
            coboundary = BarrattEccles.adem_coboundary(q, p)
            self.assertEqual(coboundary.arity, 4)
            self.assertEqual(coboundary.degree, p + q + 1)
            self.assertEqual(
                coboundary.table_reduction().boundary(),
                coboundary.boundary().table_reduction(),
            )
            residual = coboundary.boundary()
            residual += symmetric_adem_term(q, p)
            residual += symmetric_adem_term(p, q)
            self.assertFalse(residual.orbit(), msg=(q, p))

        with self.assertRaises(TypeError):
            BarrattEccles.adem_coboundary(1.0, 0)
        with self.assertRaises(TypeError):
            BarrattEccles.adem_coboundary(True, 0)
        with self.assertRaises(ValueError):
            BarrattEccles.adem_coboundary(0, 0)
        with self.assertRaises(ValueError):
            BarrattEccles.adem_coboundary(1, -1)

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
