import unittest

from comch.barratt_eccles import BarrattEcclesElement
from comch.steenrod import Steenrod
from comch.surjection import SurjectionElement
from comch.symmetric import SymmetricRing


class TestSteenrod(unittest.TestCase):
    def test_surjection_may_structure(self):
        arity = 6
        t = SymmetricRing.transposition_element(arity)
        x = Steenrod.may_structure(
            arity, 3, model='surjection').boundary()
        y = t * Steenrod.may_structure(arity, 2, model='surjection')
        self.assertEqual(x, y)

        arity = 5
        n = SymmetricRing.norm_element(arity, torsion=7)
        x = Steenrod.may_structure(
            arity, 4, model='surjection', torsion=7).boundary()
        y = n * Steenrod.may_structure(
            arity, 3, model='surjection', torsion=7)
        self.assertEqual(x, y)

        t = SymmetricRing.transposition_element(arity)
        x = Steenrod.may_structure(
            arity, 3, model='surjection',
            convention='McClure-Smith').boundary()
        y = t * Steenrod.may_structure(
            arity, 2, model='surjection',
            convention='McClure-Smith')
        self.assertEqual(x, y)

    def test_barratt_eccles_may_structure(self):
        t = SymmetricRing.transposition_element(6)
        x = Steenrod.may_structure(
            6, 3, model='barratt_eccles').boundary()
        y = t * Steenrod.may_structure(6, 2, model='barratt_eccles')
        self.assertEqual(x, y)

        n = SymmetricRing.norm_element(5, torsion=7)
        x = Steenrod.may_structure(
            5, 6, model='barratt_eccles', torsion=7).boundary()
        y = n * Steenrod.may_structure(
            5, 5, model='barratt_eccles', torsion=7)
        self.assertEqual(x, y)

    def test_surjection_operation(self):
        x = Steenrod.operation(3, 2, 3)
        y = SurjectionElement({
            (1, 2, 3, 1, 2): 2,
            (1, 3, 1, 2, 3): 2,
            (1, 2, 3, 2, 3): 2,
        })
        self.assertEqual(x, y)

    def test_barratt_eccles_operation(self):
        x = Steenrod.operation(3, 2, 3, model='barratt_eccles')
        y = BarrattEcclesElement({
            ((1, 2, 3), (2, 3, 1), (3, 1, 2)): 2,
            ((1, 2, 3), (3, 1, 2), (1, 2, 3)): 2,
        })
        self.assertEqual(x, y)

    def test_inputs(self):
        with self.assertRaises(ValueError):
            Steenrod.operation(2, 0, 0, model='not-a-model')
        with self.assertRaises(ValueError):
            Steenrod.operation(
                2, 0, 0, model='barratt_eccles',
                convention='McClure-Smith')
        with self.assertRaises(ValueError):
            Steenrod.chain(2, -1, -2, shape='not-a-shape')


if __name__ == '__main__':
    unittest.main()
