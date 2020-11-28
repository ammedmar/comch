from module import Module_element
from itertools import product, combinations


class Cube(tuple):
    '''...

    '''
    intervals: tuple

    def __init__(self, data):

        self.intervals = tuple(idx for idx, x in enumerate(data) if x == 2)

        super(Cube, self).__init__()

    @property
    def dimension(self):
        return self.count(2)

    def face(self, i, epsilon):
        '''...'''
        idx = self.intervals[i]
        answer = self[:idx] + ((epsilon + 1) % 2, ) + self[idx + 1:]
        return Cube(answer)


class CubicalEilenbergZilber_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        '''...'''

        if data:
            new_data = {}
            for k, v in data.items():
                new_k = tuple(Cube(cube) for cube in k)
                new_data[new_k] = v
            data = new_data

        super(CubicalEilenbergZilber_element, self).__init__(
            data=data, torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',').replace('2', 'e')

    @property
    def arity(self):
        '''...'''
        arities = set(len(multicube) for multicube in self.keys())
        if len(arities) != 1:
            return None
        return arities.pop()

    @property
    def degree(self):
        '''...'''
        degs = {sum(cube.dimension for cube in k) for k in self.keys()}
        if len(degs) != 1:
            return None
        return degs.pop()

    def boundary(self):
        '''Boundary of an element in a tensor product of the standard
        chains.

        # squares to zero

        >>> elmt = CubicalEilenbergZilber_element({((0, 2), (2, 1), (2, 0)): 1})
        >>> print(elmt.boundary().boundary())
        0

        '''
        answer = self.zero()
        for k, v in self.items():
            for idx, cube in enumerate(k):
                acc_dim = sum((cube_l.dimension for cube_l in k[:idx]))
                for i in range(cube.dimension):
                    for epsilon in (0, 1):
                        new_cube = cube.face(i, epsilon)
                        new_k = k[:idx] + (new_cube,) + k[idx + 1:]
                        sign_exp = (acc_dim + i + epsilon) % 2
                        answer += answer.create({new_k: v * (-1)**sign_exp})
        return answer

    def iterated_diagonal(self, n=1):
        '''Serre chain approximation to the diagonal applied n-times.

        Examples
        --------

        # chain map check:

        >>> x = CubicalEilenbergZilber_element({((0, 2, 2, 1, 2),): 1})
        >>> d_delta_x = x.iterated_diagonal(2).boundary()
        >>> delta_d_x = x.boundary().iterated_diagonal(2)
        >>> d_delta_x == delta_d_x
        True

        '''
        def sign(p):
            '''Counts the number of pairs appearing in reversed order.'''
            to_count = filter(lambda x: x[0] > x[1], combinations(p, 2))
            sign_exp = sum(1 for _ in to_count) % 2
            return (-1)**sign_exp

        def splitted(fixed, i):
            '''...'''
            if i < fixed:
                return 0
            elif i == fixed:
                return 2
            else:
                return 1

        if self.degree is None:
            raise TypeError(f'only for homogeneous elements')

        if self.arity != 1:
            raise TypeError(f'only for arity 1 elements')

        answer = self.zero()
        for k, v in self.items():
            cube = k[0]
            intervals = cube.intervals
            base = [i for idx, i in enumerate(cube) if idx not in intervals]
            for p in product(range(n), repeat=self.degree):
                multibase = [list(base) for _ in range(n)]
                for idx, fixed in enumerate(p):
                    at = intervals[idx]
                    for i, new_base in enumerate(multibase):
                        to_insert = splitted(fixed, i)
                        new_base.insert(at, to_insert)
                new_k = tuple(Cube(x) for x in multibase)
                answer += answer.create({new_k: v * sign(p)})
        return answer


if __name__ == "__main__":
    import doctest
    doctest.testmod()
