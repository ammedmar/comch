from ..basics import Module_element, TorsionError
from ..basics import SymmetricGroup_element, ArityError
from ..basics import SymmetricRing_element, SymmetricRing

from ..surjection import Surjection_element
from ..utils import partitions
from itertools import chain, product
from math import factorial


class BarrattEccles_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):

        # check input data: dict with tuple of tuple of int keys
        if data:
            if not (isinstance(data, dict)
                    and all(isinstance(x, tuple) for x in data.keys())
                    and all(isinstance(perm, tuple) for perm in
                            chain.from_iterable(data.keys()))
                    and all(isinstance(i, int) for i in
                            chain.from_iterable(chain.from_iterable(data.keys())))
                    ):
                raise TypeError('data type must be dict '
                                + 'with tuple of tuple of int keys')

            if any((set(perm) != set(range(1, len(perm) + 1)) for perm in
                    chain.from_iterable(data.keys()))):
                raise TypeError('keys must tuples of '
                                + 'permutations of (1,2,...,r)')

            # transform tuples to symmetric group elements
            new_data = {}
            for k, v in data.items():
                new_key = tuple(SymmetricGroup_element(pi) for pi in k)
                new_data[new_key] = v
        else:
            new_data = data

        # initializing element
        super(BarrattEccles_element, self).__init__(data=new_data,
                                                    torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    @property
    def arity(self):
        '''Arity of the Barratt-Eccles element, set to None if non-homogeneous
        or zero

        >>> b = BarrattEccles_element({((1, 3, 2), (2, 3, 1)): 1})
        >>> print(b.arity)
        3

        '''
        if not self:
            return None
        arities = set()
        for k in self.keys():
            arities_in_k = set(max(perm) for perm in k)
            arities |= arities_in_k  # in place union
        if len(arities) > 1:
            return None
        return arities.pop()

    @property
    def degree(self):
        pass

    @property
    def complexity(self):
        pass

    def boundary(self):
        '''Boundary as normalized chains of a simplicial set.

        >>> b = BarrattEccles_element({((1, 3, 2), (2, 3, 1)): 1})
        >>> print(b.boundary().boundary())
        0

        '''

        sign = {0: 1, 1: -1}
        answer = self.zero()
        for spx, coeff in self.items():
            for i in range(len(spx)):
                i_term = {tuple(spx[: i] + spx[i + 1:]): sign[i % 2] * coeff}
                to_add = answer.create(i_term)
                answer += to_add
        answer._reduce_rep()

        return answer

    def _reduce_rep(self):
        '''deletes degenerate keys.

        '''
        for simplex, v in self.items():
            for i in range(len(simplex) - 1):
                if simplex[i] == simplex[i + 1]:
                    self[simplex] = 0

        super()._reduce_rep()

    def __rmul__(self, other):
        '''Left action by the appropriate symmetric group ring.

        # >>> surj = Surjection_element({(1, 2, 3, 1, 2): 1})
        # >>> print(-1 * surj)
        # - (1,2,3,1,2)
        '''

        if isinstance(other, int):
            return super().__rmul__(other)

        if not isinstance(other, SymmetricRing_element):
            raise NotImplementedError

        if self.torsion != other.torsion:
            raise TorsionError

        if self.arity != other.arity:
            raise ArityError

        answer = self.zero()
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            new_k = tuple(k2 * pi for pi in k1)
            answer += self.create({new_k: v1 * v2})
        return answer

    def compose(self, *others):
        '''WRONG, NOT DONE BY ME. NOT A CHAIN MAP'''

        def _paths(p, q):
            '''returns as a list all increasing paths from (0,0) to (p,q).'''

            if (p, q) == (0, 0):
                return [((0, 0),)]
            answer = list()
            if p > 0:
                west = _paths(p - 1, q)
                for path in west:
                    answer.append(path + ((p, q),))
            if q > 0:
                south = _paths(p, q - 1)
                for path in south:
                    answer.append(path + ((p, q),))
            return answer

        def _sgn_of_path(path):
            '''...'''
            segments = range(1, len(path))
            horizontal_segments, vertical_segments = [], []
            for i in segments:
                vertex1, vertex2 = path[i - 1], path[i]
                if vertex2[0] > vertex1[0]:
                    horizontal_segments.append(i)
                else:
                    vertical_segments.append(i)
            ordered_segments = horizontal_segments + vertical_segments
            # find the permutation that transforms segments to orderedSegments
            permutation = {}
            for seg in segments:
                for j in range(1, len(ordered_segments) + 1):
                    if seg == ordered_segments[j - 1]:
                        permutation[seg] = j
            # compute the sign of the permutation
            sgn = 1
            for i in range(1, len(segments) + 1):
                for j in range(i + 1, len(segments) + 1):
                    diff = permutation[j] - permutation[i]
                    sgn *= diff // abs(diff)
            return sgn

        # partial composition
        if len(others) == 2 and isinstance(others[1], int):
            other, k = others
            if self.torsion != other.torsion:
                raise TypeError('not the same torsion')

            answer = self.zero()
            for perm_vect1, coeff1 in self.items():
                for perm_vect2, coeff2 in other.items():
                    comp = BarrattEccles_element().copy_attrs_from(answer)
                    p, q = len(perm_vect1) - 1, len(perm_vect2) - 1
                    # summands parametrized by paths from (0,0) to (p,q)
                    for path in _paths(p, q):
                        new_perm_vect = ()
                        for i, j in path:
                            perm1 = SymmetricRing_element(
                                {perm_vect1[i]: 1}, torsion=self.torsion)
                            perm2 = SymmetricRing_element(
                                {perm_vect2[j]: 1}, torsion=self.torsion)
                            partial_comp = perm1.compose(perm2, k)
                            new_perm_vect += (tuple(partial_comp.keys())[0],)
                        sgn = _sgn_of_path(path)
                        comp += BarrattEccles_element(
                            {new_perm_vect: sgn}, torsion=self.torsion)
                    answer += coeff1 * coeff2 * comp
            return answer

        # total composition
        else:
            if not len(others) == self.arity:
                raise TypeError('the number of arguments must be equal to '
                                + 'the arity of self')
            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer.compose(other, idx + 1)
            return answer

    def table_reduction(self):
        '''Returns the image of the table refuction morphism applied to self.

        # Berger-Fresse Example before Theorem 1.3.2:

        >>> b = BarrattEccles_element({((1,2,3,4), (1,4,3,2), (1,2,4,3)): 1})
        >>> print(b.table_reduction())
        (1,2,4,2,4,3) + (1,2,4,3,2,3)

        # Chain map check:

        >>> b = BarrattEccles_element({((1,2,3,4), (1,4,3,2)): 1, \
                                       ((1,2,4,3), (3,4,2,1)): 2})
        >>> dTR_b = b.table_reduction().boundary()
        >>> TRd_b = b.boundary().table_reduction()
        >>> dTR_b == TRd_b
        True

        '''

        answer = Surjection_element(torsion=self.torsion,
                                    convention='Berger-Fresse')
        for k1, v in self.items():
            d, a = len(k1) - 1, max(k1[0])
            for pi in partitions(d + a, d + 1, ordered=True):
                k2, removed = [], []
                degenerate = False
                for idx, i in enumerate(pi):
                    filtered = [i for i in k1[idx]
                                if i not in removed]
                    if idx > 0 and k2[-1] == filtered[0]:
                        degenerate = True
                        break
                    if i > 1:
                        removed += filtered[: i - 1]
                    k2 += filtered[: i]

                if not degenerate:
                    answer += answer.create({tuple(k2): v})

        answer._reduce_rep()

        return answer

    def orbit(self, representation='trivial'):
        '''Orbit under the symmetric group action

        SIGN PERMUTATIONS FAILING

        '''
        if not self:
            return self

        answer = BarrattEccles_element(torsion=self.torsion)
        for k, v in self.items():
            inverse = tuple(k[0].index(i + 1) + 1 for i in range(len(k[0])))
            permutation = SymmetricRing_element({inverse: 1}, torsion=self.torsion)
            if representation == 'sign':
                permutation = sign(k[0]) * permutation
            answer += permutation * BarrattEccles_element({k: v}, torsion=self.torsion)

        return answer

    def alexander_whitney(self, r=1):
        '''...'''

        def split(multispx):
            a, b = multispx[0], multispx[1:]
            return set((a[:i + 1], a[i:]) + b for i in range(len(a)))

        answer = Module_element(torsion=self.torsion)
        for k, v in self.items():
            to_add = set(((k,),))
            for s in range(1, r + 1):
                to_add = set.union(*(split(multispx) for multispx in to_add))
            answer += Module_element(
                {multispx: v for multispx in to_add}).copy_attrs_from(self)

        return answer


class BarrattEccles():
    '''Class producing Barratt-Eccles elements of special interest.'''

    @staticmethod
    def steenrod_product(arity, degree, torsion=None):
        '''Returns a Barratt-Eccles element representing the Steenrod
        product in the given arity and degree.

        # chain map checks:

        >>> t = SymmetricRing.transposition_element(6)
        >>> x = BarrattEccles.steenrod_product(6, 3).boundary()
        >>> y = t * BarrattEccles.steenrod_product(6, 2)
        >>> print(x == y)
        True

        >>> n = SymmetricRing.norm_element(5, torsion=7)
        >>> x = BarrattEccles.steenrod_product(5, 6, torsion=7).boundary()
        >>> y = n * BarrattEccles.steenrod_product(5, 5, torsion=7)
        >>> print(x == y)
        True

        '''

        operators = {
            0: SymmetricRing.norm_element(arity),
            1: SymmetricRing.transposition_element(arity)
        }

        def _psi(arity, degree):
            '''Recursive definition over the integers.'''

            if degree == 0:
                return BarrattEccles_element(
                    {(tuple(range(1, arity + 1)),): 1})
            else:
                previous = _psi(arity, degree - 1)
                acted_on = operators[degree % 2] * previous
                identity = tuple(range(1, arity + 1))
                new_data = {(identity,) + k: v for k, v in acted_on.items()}
                return BarrattEccles_element(new_data)

        integral_answer = _psi(arity, degree)
        if torsion:
            integral_answer.set_torsion(torsion)
        return integral_answer
