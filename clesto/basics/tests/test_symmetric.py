import unittest
from clesto.basics import SymmetricGroup_element, SymmetricRing_element


class TestSymmetricGroup(unittest.TestCase):

    def setUp(self):
        self.x = SymmetricGroup_element((5, 2, 4, 3, 1))
        self.y = SymmetricGroup_element((1, 2, 3, 4, 5))

    def test_sign(self):
        self.assertEqual(self.x.sign, 1)

    def test_arity(self):
        self.assertEqual(self.x.arity, 5)

    def test_inverse(self):
        self.assertEqual(self.x.inverse(), self.x)

    def test_mul(self):
        self.assertEqual(self.x * self.y, self.x)

    def test_pow(self):
        self.assertEqual(self.x ** 2, self.y)

    def test_compose(self):
        self.assertEqual(self.x.compose(SymmetricGroup_element((2, 1)), 3),
                         SymmetricGroup_element((6, 2, 5, 4, 3, 1)))


class TestSymmetricRing(unittest.TestCase):

    def setUp(self):
        self.x = SymmetricRing_element({(5, 2, 4, 3, 1): 1})
        self.y = SymmetricRing_element({(1, 2, 3, 4, 5): 1})

    def test_arity(self):
        self.assertEqual(self.x.arity, 5)

    def test_mul(self):
        self.assertEqual(3 * self.x,
                         SymmetricRing_element({(5, 2, 4, 3, 1): 3}))

        self.assertEqual(self.x * self.y, self.x)

    def test_pow(self):
        self.assertEqual(self.x ** 2, self.y)

    def test_compose(self):
        self.assertEqual(self.x.compose(SymmetricRing_element({(2, 1): 1}), 3),
                         SymmetricRing_element({(6, 2, 5, 4, 3, 1): 1}))


if __name__ == '__main__':
    unittest.main()
