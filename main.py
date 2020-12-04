# Run tests
# ---------

import doctest
import clesto
doctest.testmod(clesto.basics.module)
doctest.testmod(clesto.basics.symmetric)
doctest.testmod(clesto.barratt_eccles.barratt_eccles)
doctest.testmod(clesto.surjection.surjection)
doctest.testmod(clesto.eilenberg_zilber.simplicial)
doctest.testmod(clesto.eilenberg_zilber.cubical)
print('tests passed')

# ------------------------------------------------------------------------------
from clesto import *
# x = Surjection_element({(1, 2, 4, 3): 1}, convention='McClure-Smith', torsion=3)
# print(x.zero().convention)
