import unittest
from comch.symmetric import SymmetricGroupElement, SymmetricRingElement


class TestSymmetricGroup(unittest.TestCase):

    def setUp(self):
        self.x = SymmetricGroupElement((5, 2, 4, 3, 1))
        self.y = SymmetricGroupElement((1, 2, 3, 4, 5))

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
        self.assertEqual(self.x.compose(SymmetricGroupElement((2, 1)), 3),
                         SymmetricGroupElement((6, 2, 5, 4, 3, 1)))


class TestSymmetricRing(unittest.TestCase):

    def setUp(self):
        self.x = SymmetricRingElement({(5, 2, 4, 3, 1): 1})
        self.y = SymmetricRingElement({(1, 2, 3, 4, 5): 1})

    def test_arity(self):
        self.assertEqual(self.x.arity, 5)

    def test_mul(self):
        self.assertEqual(3 * self.x,
                         SymmetricRingElement({(5, 2, 4, 3, 1): 3}))

        self.assertEqual(self.x * self.y, self.x)

    def test_pow(self):
        self.assertEqual(self.x ** 2, self.y)

    def test_compose(self):
        self.assertEqual(self.x.compose(SymmetricRingElement({(2, 1): 1}), 3),
                         SymmetricRingElement({(6, 2, 5, 4, 3, 1): 1}))


if __name__ == '__main__':
    unittest.main()
