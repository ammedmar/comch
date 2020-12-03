from ..basics import Module_element, TorsionError
from ..basics import SymmetricModule_element, ArityError
from ..utils import pairwise
from itertools import chain, product, combinations_with_replacement


class Simplex(tuple):
    '''...'''

    @property
    def dimension(self):
        '''...'''

        return len(self) - 1

    def face(self, i):
        '''...'''

        return Simplex(self[:i] + self[i + 1:])

    def degeneracy(self, i):
        '''...'''

        return Simplex(self[:i + 1] + self[i:])

    def coface(self, i):
        '''...'''

        def d_i(i, j): return j + 1 if j >= i else j

        return tuple(d_i(i, j) for j in self)

    def codegeneracy(self, i):
        '''...'''

        def s_i(i, j): return j - 1 if j > i else j

        return tuple(s_i(i, j) for j in self)

    def is_degenerate(self):
        '''...'''

        conseq_values = any([i == j for i, j in pairwise(self)])
        empty_simplex = (self.dimension == -1)
        return empty_simplex or conseq_values


class EilenbergZilber_element(Module_element):
    '''...'''

    def __init__(self, data=None, dimension=None, torsion=None):
        '''...'''

        if data:
            if not dimension:
                dimension = max(i for i in chain.from_iterable(
                                chain.from_iterable(data.keys())))

            new_data = {}
            for k, v in data.items():
                new_k = tuple(Simplex(spx) for spx in k)
                new_data[new_k] = v
            data = new_data

        self.dimension = dimension

        super(EilenbergZilber_element, self).__init__(data=data,
                                                      torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def _latex_(self):
        '''Representation in Latex.

        Example
        -------

        >>> x = EilenbergZilber_element({((0,), (0, 1, 2)): 1, \
                                       ((0, 1), (1, 2)): 1, \
                                       ((0, 1, 2), (2,)): 1})
        >>> print(x._latex_())
        [0] \otimes [0,1,2] + [0,1] \otimes [1,2] + [0,1,2] \otimes [2]

        '''
        string = str(self)
        string = string.replace(',),(', '] \otimes [')
        string = string.replace('),(', '] \otimes [')
        string = string.replace('((', '[')
        string = string.replace(',))', ']').replace('),)', ']')
        string = string.replace('))', ']')

        return string

    def create(self, data=None):
        '''...'''
        return type(self)(data, torsion=self.torsion,
                          dimension=self.dimension)

    def zero(self):
        '''...'''
        return self.create()

    @property
    def arity(self):
        arities = set(len(multispx) for multispx in self.keys())
        if len(arities) != 1:
            return None
        return arities.pop()

    @property
    def degree(self):
        '''...'''
        degs = {sum(spx.dimension for spx in k) for k in self.keys()}
        if len(degs) != 1:
            return None
        return degs.pop()

    def boundary(self):
        '''Boundary of an element in a tensor product of the standard
        chains.

        # squares to zero

        >>> elmt = EilenbergZilber_element({((0, 1), (0,), (0, 1, 2, 3)): 1})
        >>> print(elmt.boundary().boundary())
        0

        '''
        answer = self.zero()
        for k, v in self.items():
            for idx, spx in enumerate(k):
                acc_dim = sum((spx_l.dimension for spx_l in k[:idx]))
                for i in range(spx.dimension + 1):
                    new_spx = spx.face(i)
                    new_k = k[:idx] + (new_spx,) + k[idx + 1:]
                    sign_exp = (acc_dim + i) % 2
                    answer += answer.create({new_k: v * (-1)**sign_exp})
        return answer

    def __rmul__(self, other):
        '''Left action by the appropriate symmetric group ring.

        # chain map checks:

        >>> rho = SymmetricModule_element({(2, 3, 1): 1})
        >>> elmt = EilenbergZilber_element({((0, 1), (0,), (0, 1, 2, 3)): 1})
        >>> x, y = (rho * elmt).boundary(), rho * elmt.boundary()
        >>> x == y
        True

        '''
        def sign(perm, multispx):
            signs = {0: 1, 1: -1}
            weights = [spx.dimension % 2 for spx in k1]
            answer = 0
            for idx, i in enumerate(perm):
                right = [weights[perm.index(j)] for
                         j in perm[idx + 1:] if i > j]
                answer += sum(right) * weights[idx]
            return signs[answer % 2]

        if isinstance(other, int):
            return super().__rmul__(other)

        if not isinstance(other, SymmetricModule_element):
            raise TypeError(f'right mult. by type int or \
                SymmetricModule_element not {type(other)}')

        if self.torsion != other.torsion:
            raise TorsionError

        if self.arity != other.arity:
            raise ArityError

        answer = self.zero()
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            new_key = tuple(k1[k2.index(i + 1) - 1] for i in range(self.arity))
            new_sign = sign(k2, k1)
            answer += self.create({new_key: new_sign * v1 * v2})
        return answer

    def codegeneracy(self, i):
        '''Covariant action of the i-th codegeneracy. It is the zero map
        on any element in the Eilenberg-Zilber operad.

        '''
        if i > self.dimension - 1:
            raise TypeError('codegeneracy out of range')

        answer = self.zero()
        for k, v in self.items():
            new_k = tuple(spx.codegeneracy(i) for spx in k)
            answer += answer.create({new_k: v})
        return answer

    def coface(self, i):
        '''Covariant action of the i-th coface.'''

        if i > self.dimension:
            raise TypeError('coface out of range')

        answer = self.zero()
        for k, v in self.items():
            new_k = tuple(tuple(spx.coface(i) for spx in k))
            answer += answer.create({new_k: v})
        return answer

    def _reduce_rep(self):
        '''Sets to 0 the summands with degenerate simplices.'''

        for k, v in self.items():
            if any([spx.is_degenerate() for spx in k]):
                self[k] = 0

        super()._reduce_rep()

    def iterated_diagonal(self, n=1):
        '''Alexander-Whitney chain approximation to the diagonal applied
        n-times.

        Examples
        --------

        # chain map check:

        >>> x = EilenbergZilber_element({((0, 1, 2), ): 1})
        >>> dx = x.boundary()
        >>> dx.iterated_diagonal(3) == x.iterated_diagonal(3).boundary()
        True

        '''

        if self.degree is None:
            raise TypeError(f'only for homogeneous elements')

        if self.arity != 1:
            raise TypeError(f'only for arity 1 elements')

        answer = self.zero()
        for k, v in self.items():
            spx = k[0]
            for p in combinations_with_replacement(range(self.degree + 1), n):
                p = (0,) + p + (self.degree,)
                new_k = []
                for i, j in pairwise(p):
                    new_k.append(Simplex(spx[i:j + 1]))
                answer += self.create({tuple(new_k): v})
        return answer


class EilenbergZilber():
    '''Class producing Eilenberg-Zilber elements of special interest.'''

    def standard_element(n, torsion=None):
        return EilenbergZilber_element({(tuple(range(n + 1)), ): 1},
                                       torsion=torsion)

    def boundary_element(n, torsion=None):
        '''...'''

        sign = {0: 1, 1: -1}
        answer = EilenbergZilber_element(dimension=n, torsion=torsion)
        spx = Simplex(range(n + 1))
        for i in range(n + 1):
            new_k = (spx[:i] + spx[i + 1:], )
            answer += answer.create({new_k: sign[i % 2]})
        return answer

    def coproduct_element(n, torsion=None):
        '''...'''

        answer = EilenbergZilber_element(dimension=n, torsion=torsion)
        spx = Simplex(range(n + 1))
        for i in range(1, n + 2):
            new_k = (spx[:i], spx[i - 1:])
            answer += answer.create({new_k: 1})
        return answer
