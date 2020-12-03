# Run tests
# ---------

# import doctest
# import clesto
# doctest.testmod(clesto.module)
# doctest.testmod(clesto.symmetric)
# doctest.testmod(clesto.barratt_eccles)
# doctest.testmod(clesto.surjection)
# doctest.testmod(clesto.simplicial)
# doctest.testmod(clesto.cubical)
# print('tests passed')

# ------------------------------------------------------------------------------
from clesto import *

x = CubicalEilenbergZilber_element({((2,), (1,)): 1, ((0,), (2,)): 1})
print(x._latex_())
