import clesto
test = False
if test:
    import doctest
    doctest.testmod(clesto.module)
    doctest.testmod(clesto.symmetric)
    doctest.testmod(clesto.barratt_eccles)
    doctest.testmod(clesto.surjection)
    doctest.testmod(clesto.simplicial)
    doctest.testmod(clesto.cubical)

    print('good')

#------------------------------------------------------------------------------
