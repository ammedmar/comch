from ..basics import Module_element
from ..utils import pairwise
from itertools import combinations, product


class Cube(tuple):
    """Models an elementary cube"""
    intervals: tuple

    def __init__(self, data):
        """Initialize"""

        self.intervals = tuple(idx for idx, x in enumerate(data) if x == 2)

        super(Cube, self).__init__()

    @property
    def dimension(self):
        return self.count(2)

    def face(self, i, epsilon):
        """..."""
        idx = self.intervals[i]
        answer = self[:idx] + ((epsilon + 1) % 2, ) + self[idx + 1:]
        return Cube(answer)


class CubicalEilenbergZilber_element(Module_element):
    """..."""

    dimenion: int = None

    def __init__(self, data=None, torsion=None):
        """..."""

        if data:
            new_data = {}
            for k, v in data.items():
                dim_k = set(len(cube) for cube in k)
                if len(dim_k) != 1:
                    raise TypeError('some faces have different length')
                new_k = tuple(Cube(cube) for cube in k)
                new_data[new_k] = v
            data = new_data
            dimension = dim_k.pop()  # they are all equal, using last.
            self.dimension = dimension

        super(CubicalEilenbergZilber_element, self).__init__(
            data=data, torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def _latex_(self):
        """Representation in Latex.

        ERROR: (2,2) --> [01]2 since no coma is left

        >>> x = CubicalEilenbergZilber_element({((2,), (1,)): 1, ((0,), (2,)): 1})
        >>> print(x._latex_())
        [01] \otimes [1] + [0] \otimes [01]

        """
        string = str(self)
        string = string.replace('1', '[1]').replace('0', '[0]')
        string = string.replace('2,', '[01]').replace(',2', '[01]')
        string = string.replace(',),(', ' \otimes ')
        string = string.replace('),(', ' \otimes ')
        string = string.replace(')', '').replace('((', '')
        string = string.replace(')', '').replace('))', '')
        string = string.replace(',', '')

        return string

    @property
    def arity(self):
        """..."""
        arities = set(len(multicube) for multicube in self.keys())
        if len(arities) != 1:
            return None
        return arities.pop()

    @property
    def degree(self):
        """..."""
        degs = {sum(cube.dimension for cube in k) for k in self.keys()}
        if len(degs) != 1:
            return None
        return degs.pop()

    def boundary(self):
        """Boundary of an element in a tensor product of the standard
        chains.

        # squares to zero

        >>> elmt = CubicalEilenbergZilber_element({((0, 2), (2, 1), (2, 0)): 1})
        >>> print(elmt.boundary().boundary())
        0

        """
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
        """Serre chain approximation to the diagonal applied n-times.

        Examples
        --------

        # chain map check:

        >>> x = CubicalEilenbergZilber_element({((0, 2, 2, 1, 2),): 1})
        >>> d_delta_x = x.iterated_diagonal(2).boundary()
        >>> delta_d_x = x.boundary().iterated_diagonal(2)
        >>> d_delta_x == delta_d_x
        True

        """
        def sign(p):
            """Counts the number of pairs appearing in reversed order.

            """
            to_count = filter(lambda x: x[0] > x[1], combinations(p, 2))
            sign_exp = sum(1 for _ in to_count) % 2
            return (-1)**sign_exp

        def elementary_summand(fixed, i):
            """Models as a function the element 0,...,0,2,1,...,1 appearing
            as one of the summands of the iterated diagonal of an interval.

            """
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
            for p in product(range(n + 1), repeat=self.degree):
                multibase = [list(base) for _ in range(n + 1)]
                for idx, fixed in enumerate(p):
                    at = intervals[idx]
                    for i, new_base in enumerate(multibase):
                        to_insert = elementary_summand(fixed, i)
                        new_base.insert(at, to_insert)
                new_k = tuple(Cube(x) for x in multibase)
                answer += answer.create({new_k: v * sign(p)})
        return answer

    def join(self):
        """Join of an element in the cubical EZ operad thought of
        as an element in the tensor product, computed using the left
        comb. (I suspect it is associative product though.)

        Examples
        --------

        # boundary of the join

        >>> x = CubicalEilenbergZilber_element({((0, 0, 1), \
                                                 (1, 0, 0)): 1})
        >>> print(x.join().boundary() + x.boundary().join())
        ((1,0,0),) - ((0,0,1),)

        """
        def is_zero(left, right):
            """Two conditions need to be satisfied for a triple
            (cube1, cube2, i) give a nonzero i-join: no intervals
            in cube1 can have indices greater or equal to i, and
            no intervals in cube2 can have indices less than or
            equal to i.

            """
            if left == tuple() or right == tuple():
                return False
            if isinstance(right, int):
                return right <= max(left)
            if isinstance(left, int):
                return left >= min(right)

        def _join(i, cube1, cube2, sign_exp):
            """the i-th elementary join keeping track of signs.

            """
            cube = Cube(cube1[:i] + (2,) + cube2[i + 1:])
            p, q = cube1[i], cube2[i]
            if (p, q) == (0, 1):
                return cube, sign_exp
            elif (p, q) == (1, 0):
                return cube, (1 + sign_exp) % 2
            else:
                return None, None

        if not self:
            return self

        if self.degree is None:
            raise TypeError(f'only for homogeneous elements')

        answer = self.zero()
        for k, v in self.items():
            for inds in combinations(range(self.dimension), self.arity - 1):
                skip = False
                for i, (cube1, cube2) in zip(inds, pairwise(k)):
                    if (is_zero(cube1.intervals, i) or
                            is_zero(i, cube2.intervals)):
                        skip = True
                        break
                if not skip:
                    non_zero = True
                    sign_exp = 0
                    cube = k[0]
                    for i, next_cube in zip(inds, k[1:]):
                        cube, sign_exp = _join(i, cube, next_cube, sign_exp)
                        if cube is None:
                            non_zero = False
                            break
                    if non_zero:
                        answer += answer.create({(cube, ): (-1)**sign_exp})

        return answer


class CubicalEilenbergZilber:
    """.."""

    def standard_element(n, torsion=None):
        """..."""
        return CubicalEilenbergZilber_element({((2,) * n, ): 1}, torsion=torsion)
