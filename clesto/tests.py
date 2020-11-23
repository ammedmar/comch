from clesto import *

# for a, b, c in product((1, -1), repeat=3):

x = EilenbergZilber_element({((0, 1, 2), (0, 1)): -1,
                             ((0, 1, 2), (1, 2)): -1,
                             ((0, 2), (0, 1, 2)): 1})

y = EilenbergZilber_element({((1, 2), (1, 2)): 1,
                             ((0, 2), (0, 2)): -1,
                             ((0, 1), (0, 1)): 1})

# print(x.boundary() + y)


z = EilenbergZilber_element({((0, 1, 2), (0,)): 1,
                             ((1, 2), (0, 1)): -1,
                             ((2, ), (0, 1, 2)): 1})

# print(z.boundary())

print(Simplex(()).dimension)

# print(x.boundary().codegeneracy(1))
# print(x.codegeneracy(1))
# print(x.coface(0))

# print(EilenbergZilber.boundary(3))
# print(EilenbergZilber.coproduct(3))