import unittest
from comch.barratt_eccles import BarrattEcclesElement
from comch.cubical import CubicalElement
from comch.free_module import FreeModuleElement
from comch.simplicial import SimplicialElement
from comch.surjection import SurjectionElement
from comch.symmetric import SymmetricRingElement


class TestModule(unittest.TestCase):

    def setUp(self):
        self.x = FreeModuleElement({'a': 3, 'b': 1})
        self.y = FreeModuleElement({'a': 2, 'b': -1})

    def test_hash(self):
        a = {self.x, self.x}
        b = {self.x}
        self.assertEqual(a, b)

    def test_add(self):
        self.assertEqual(self.x + self.y, FreeModuleElement({'a': 5}))

    def test_sub(self):
        self.assertEqual(self.x - self.y, FreeModuleElement({'a': 1, 'b': 2}))

    def test_rmul(self):
        self.c = 3
        self.assertEqual(3 * self.x, FreeModuleElement({'a': 9, 'b': 3}))

    def test_truediv(self):
        self.c = 3
        z = FreeModuleElement({'a': 1, 'b': 2}, torsion=5)
        self.assertEqual(z / 3, FreeModuleElement({'a': 2, 'b': 4}, torsion=5))

    def test_neg(self):
        self.assertEqual(- self.x, FreeModuleElement({'a': -3, 'b': -1}))

    def test_iadd(self):
        self.x += self.y
        self.assertEqual(self.x, FreeModuleElement({'a': 5}))

    def test_isub(self):
        self.x -= self.y
        self.assertEqual(self.x, FreeModuleElement({'a': 1, 'b': 2}))

    def test_torsion(self):
        class ElementWithDefaultTorsion(FreeModuleElement):
            default_torsion = -5

        self.assertEqual(
            FreeModuleElement({'a': 1}, torsion=0),
            FreeModuleElement({'a': 1}),
        )
        self.assertEqual(ElementWithDefaultTorsion().torsion, 5)
        self.assertEqual(FreeModuleElement({'a': 1}, torsion=1),
                         FreeModuleElement(torsion=1))
        self.assertEqual(
            FreeModuleElement({'a': -1}, torsion=-5),
            FreeModuleElement({'a': 4}, torsion=5),
        )
        self.assertEqual(FreeModuleElement(torsion=-5).torsion, 5)

        x = FreeModuleElement({'a': 1})
        self.assertEqual(x.set_torsion(1), FreeModuleElement(torsion=1))
        self.assertEqual(FreeModuleElement({'a': -1}).set_torsion(-5),
                         FreeModuleElement({'a': 4}, torsion=5))

        for torsion in ('free', 1.5, True):
            with self.subTest(torsion=torsion):
                with self.assertRaises(TypeError):
                    FreeModuleElement({'a': 1}, torsion=torsion)
                with self.assertRaises(TypeError):
                    FreeModuleElement({'a': 1}).set_torsion(torsion)

    def test_torsion_one_subclasses(self):
        for torsion in (-1, 1):
            elements = [
                SymmetricRingElement({(1, 2): 1}, torsion=torsion),
                BarrattEcclesElement({((1, 2),): 1}, torsion=torsion),
                SurjectionElement({(1, 2): 1}, torsion=torsion),
                SimplicialElement({((0, 1),): 1}, torsion=torsion),
                CubicalElement({((2,),): 1}, torsion=torsion),
            ]
            self.assertTrue(all(not element for element in elements))


if __name__ == '__main__':
    unittest.main()
