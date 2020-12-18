import unittest
from clesto.module import Module_element


class TestModule(unittest.TestCase):

    def setUp(self):
        self.x = Module_element({'a': 3, 'b': 1})
        self.y = Module_element({'a': 2, 'b': -1})

    def test_hash(self):
        a = {self.x, self.x}
        b = {self.x}
        self.assertEqual(a, b)

    def test_add(self):
        self.assertEqual(self.x + self.y, Module_element({'a': 5}))

    def test_sub(self):
        self.assertEqual(self.x - self.y, Module_element({'a': 1, 'b': 2}))

    def test_rmul(self):
        self.c = 3
        self.assertEqual(3 * self.x, Module_element({'a': 9, 'b': 3}))

    def test_neg(self):
        self.assertEqual(- self.x, Module_element({'a': -3, 'b': -1}))

    def test_iadd(self):
        self.x += self.y
        self.assertEqual(self.x, Module_element({'a': 5}))

    def test_isub(self):
        self.x -= self.y
        self.assertEqual(self.x, Module_element({'a': 1, 'b': 2}))


if __name__ == '__main__':
    unittest.main()
