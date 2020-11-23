from clesto.module_element import Module_element, TorsionError
# from clesto.symmetric import SymmetricGroup_element, \
#     SymmetricModule_element, SymmetricModule, ArityError
from clesto.utils import pairwise, decompositions, distinct_permutations
from itertools import chain, combinations, product
from operator import attrgetter


class Simplex(tuple):
    '''...'''

    @property
    def dimension(self):
        '''...'''

        return len(self) - 1

    def face(self, i):
        '''...'''

        return Simplex(self[:i] + self[i + 1:])

    def coface(self, i):
        '''...'''

        def d_i(i, j): return j + 1 if j >= i else j

        return tuple(d_i(i, j) for j in self)

    def codegeneracy(self, i):
        '''...'''

        def s_i(i, j): return j - 1 if j > i else j

        return tuple(s_i(i, j) for j in self)

    def is_degenerate(spx):
        '''...'''

        return any([i == j for i, j in pairwise(spx)])


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

    def create(self, data=None):
        '''...'''
        return type(self)(data, torsion=self.torsion,
                          dimension=self.dimension)

    def zero(self):
        '''...'''
        return self.create()

    @ property
    def arity(self):
        arities = set(len(multispx) for multispx in self.keys())
        if len(arities) != 1:
            return None
        return arities.pop()

    def boundary(self):
        '''...'''
        sign = {0: 1, 1: -1}
        answer = self.zero()
        for k, v in self.items():
            for idx, spx in enumerate(k):
                acc_dim = sum((spx_l.dimension for spx_l in k[:idx]))
                for i in range(spx.dimension + 1):
                    new_spx = spx.face(i)
                    new_k = k[:idx] + (new_spx,) + k[idx + 1:]
                    new_v = v * sign[acc_dim % 2] * sign[i % 2]
                    answer += answer.create({new_k: new_v})
        return answer

    def codegeneracy(self, i):
        '''...'''
        if i > self.dimension - 1:
            raise TypeError('codegeneracy out of range')

        answer = self.zero()
        for k, v in self.items():
            new_k = tuple(spx.codegeneracy(i) for spx in k)
            answer += answer.create({new_k: v})
        return answer

    def coface(self, i):
        '''...'''

        if i > self.dimension:
            raise TypeError('coface out of range')

        answer = self.zero()
        for k, v in self.items():
            new_k = tuple(tuple(spx.coface(i) for spx in k))
            answer += answer.create({new_k: v})
        return answer

    def _reduce_rep(self):
        '''...'''

        for k, v in self.items():
            if any([spx.is_degenerate() for spx in k]):
                self[k] = 0

        super()._reduce_rep()


class EilenbergZilber():
    '''Class producing Eilenberg-Zilber elements of special interest.'''

    def boundary(n, torsion=None):
        '''...'''

        sign = {0: 1, 1: -1}
        answer = EilenbergZilber_element(dimension=n, torsion=torsion)
        spx = Simplex(range(n + 1))
        for i in range(n + 1):
            new_k = (spx[:i] + spx[i + 1:], )
            answer += answer.create({new_k: sign[i % 2]})
        return answer

    def coproduct(n, torsion=None):
        '''...'''

        answer = EilenbergZilber_element(dimension=n, torsion=torsion)
        spx = Simplex(range(n + 1))
        for i in range(1, n + 2):
            new_k = (spx[:i], spx[i - 1:])
            answer += answer.create({new_k: 1})
        return answer
