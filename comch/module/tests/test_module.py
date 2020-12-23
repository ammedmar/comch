import unittest
from comch.module import ModuleElement


class TestModule(unittest.TestCase):

    def setUp(self):
        self.x = ModuleElement({'a': 3, 'b': 1})
        self.y = ModuleElement({'a': 2, 'b': -1})

    def test_hash(self):
        a = {self.x, self.x}
        b = {self.x}
        self.assertEqual(a, b)

    def test_add(self):
        self.assertEqual(self.x + self.y, ModuleElement({'a': 5}))

    def test_sub(self):
        self.assertEqual(self.x - self.y, ModuleElement({'a': 1, 'b': 2}))

    def test_rmul(self):
        self.c = 3
        self.assertEqual(3 * self.x, ModuleElement({'a': 9, 'b': 3}))

    def test_neg(self):
        self.assertEqual(- self.x, ModuleElement({'a': -3, 'b': -1}))

    def test_iadd(self):
        self.x += self.y
        self.assertEqual(self.x, ModuleElement({'a': 5}))

    def test_isub(self):
        self.x -= self.y
        self.assertEqual(self.x, ModuleElement({'a': 1, 'b': 2}))


if __name__ == '__main__':
    unittest.main()
