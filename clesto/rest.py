from clesto.module_element import Module_element
from clesto.symmetric import SymmetricGroup_element, SymmetricModule_element
from clesto.surjection import Surjection_element
from clesto.utils import partitions

from itertools import chain, combinations, product
from functools import reduce
from math import floor, factorial


class CyclicModule_element(Module_element):
    '''Modeling elements in Z/mZ[C_n]

    '''

    default_order = 'infinite'

    def __init__(self, data=None, torsion=None, order=None):

        # check input data: dict with int keys
        if data:
            if not (isinstance(data, dict) and
                    all((type(k) is int for k in data.keys()))
                    ):
                raise TypeError('data type must be dict with int keys')

        # checking input order: positive int or 'infinite'
        if order is not None:
            if not (isinstance(order, int) and order > 0 or
                    order != 'infinite'
                    ):
                raise TypeError("order must be a positive int or 'infinite'")

        # setting order
        n = order if order else type(self).default_order
        setattr(self, 'order', n)

        # initializing element
        super(CyclicModule_element, self).__init__(data=data, torsion=torsion)

    def __str__(self):
        '''...'''
        if not self:
            return '0'
        else:
            answer = ''
            for exponent, coefficient in self.items():
                if coefficient != 1 or exponent == 0:
                    answer += f'+{coefficient}q^{exponent}'
                else:
                    answer += f'+q^{exponent}'
            if answer[0] == '+':
                answer = answer[1:]

            return answer.replace('q^0', '').replace('q^1', 'q')

    def __mul__(self, other):
        '''...'''
        self.compare_attributes(other)
        answer = type(self)().copy_attrs_from(self)
        for k1, v1 in self.items():
            for k2, v2 in other.items():
                answer[k1 + k2] += v1 * v2
        answer._reduce_rep()
        return answer

    def __call__(self, other):
        '''...'''
        if isinstance(other, CyclicModule_element):
            return self.__mul__(other)

        if isinstance(other, CyclicDGModule_element):
            self.compare_attributes(other)
            answer = type(other)()
            for attr, value in other.__dict__.items():
                setattr(answer, attr, value)
            for k, v1 in self.items():
                for x, v2 in other.items():
                    y = tuple(k + i for i in x)
                    answer[y] += v1 * v2
            answer._reduce_rep()
            return answer

    def _reduce_rep(self):
        '''in place mod p reduction of the keys'''

        # reducing keys mod order
        if self.order != 'infinite':
            aux = list(self.items())
            self.clear()
            for k, v in aux:
                self[k % self.order] += v
        super()._reduce_rep()

    def set_order(self, n):
        setattr(self, 'order', n)
        self._reduce_rep()
        return(self)

    @staticmethod
    def transposition_element(torsion=None, order=None):
        '''...'''
        return CyclicModule_element({1: 1, 0: -1},
                                    torsion=torsion, order=order)

    @classmethod
    def norm_element(self, torsion=None, order=None):
        '''...'''
        if order is None:
            order = self.default_order

        if order == 'infinite':
            raise TypeError('norm element not define for infinite order')

        return CyclicModule_element({i: 1 for i in range(order)},
                                    torsion=torsion, order=order)

    def psi(self, d):
        '''...'''
        # recursive function to find psi(e_d)
        def _psi_on_generator(d):
            if d <= 0:
                return CyclicDGModule_element({(0,): 1}, torsion=self.torsion,
                                              order=self.order)
            else:
                operators = {
                    0: CyclicModule_element.norm_element(torsion=self.torsion,
                                                         order=self.order),
                    1: CyclicModule_element.transposition_element(
                        torsion=self.torsion, order=self.order)}

                op_psi = operators[d % 2](_psi_on_generator(d - 1))

                return CyclicDGModule_element(
                    {(0,) + k: v for k, v in op_psi.items()},
                    torsion=self.torsion, order=self.order)

        # Compute psi knowing psi(e_d)
        answer = CyclicDGModule_element().copy_attrs_from(self)
        for k1 in _psi_on_generator(d).keys():
            for k2, v2 in self.items():
                to_add = CyclicDGModule_element(
                    {tuple(k2 + i for i in k1): v2},
                    torsion=self.torsion, order=self.order)
                answer += to_add
        return answer

    def as_BarrattEccles_element(self, d):
        '''...'''
        return self.psi(d)

    def as_Surjection_element(self, d):
        '''Uses the direct formula not the one involving the
        table_reduction morphism'''

        m = self.torsion
        assert isinstance(m, int), 'requires finite torsion'

        n = int(d / 2)
        start = tuple((0, 1))[: (d % 2) + 1]
        shifted_keys = []
        for prod in product(range(0, m), repeat=n):
            middle = tuple()
            for i in prod:
                middle += (i, (i + 1) % m)
            key = start + middle
            key += tuple((i + key[-1]) % m for i in range(1, m))
            shifted_keys.append(key)

        keys = [tuple(i + 1 for i in k) for k in shifted_keys]

        rhos = (v * SymmetricModule_element.rho(m, k) for k, v in self.items())
        # sum of all rhos
        coeff = reduce(lambda x, y: x + y, rhos).set_torsion(m)

        return coeff * Surjection_element({key: 1 for key in keys}, torsion=m)


class DGModule_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):

        # check input data: dict with tuple keys
        if data:
            if not all((isinstance(x, tuple) for x in data.keys())):
                raise ValueError('data type must be dict with tuple keys')

        # initializing element
        super(DGModule_element, self).__init__(data=data, torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def _reduce_rep(self):
        '''deletes degenerate keys and reduces as Module_element'''

        # removes degenerate simplices
        for simplex, v in self.items():
            for i in range(len(simplex) - 1):
                if simplex[i] == simplex[i + 1]:
                    self[simplex] = 0

        super()._reduce_rep()

    def boundary(self):
        '''...'''
        sign = {0: 1, 1: -1}
        bdry = type(self)().copy_attrs_from(self)
        for spx, coeff in self.items():
            for i in range(len(spx)):
                i_term = {tuple(spx[: i] + spx[i + 1:]): sign[i % 2] * coeff}
                to_add = type(self)(i_term).copy_attrs_from(bdry)
                bdry += to_add
        bdry._reduce_rep()

        return bdry

    def alexander_whitney(self, r=1):
        '''...'''

        def split(multispx):
            a, b = multispx[0], multispx[1:]
            return set((a[:i + 1], a[i:]) + b for i in range(len(a)))

        answer = Module_element().copy_attrs_from(self)
        for k, v in self.items():
            to_add = set(((k,),))
            for s in range(1, r + 1):
                to_add = set.union(*(split(multispx) for multispx in to_add))
            answer += Module_element(
                {multispx: v for multispx in to_add}).copy_attrs_from(self)

        return answer


class CyclicDGModule_element(DGModule_element):
    '''...'''

    default_order = 'infinite'

    def __init__(self, data=None, torsion=None, order=None):

        # check input data: dict with tuple of int keys
        if data:
            if not (isinstance(data, dict)
                    and all((isinstance(i, int) for i in
                             chain.from_iterable(data.keys())))
                    ):
                raise ValueError('data type must be dict'
                                 + 'with tuple of int keys')
        # set order
        n = order if order else type(self).default_order
        setattr(self, 'order', n)

        # initialize element
        super(CyclicDGModule_element, self).__init__(data=data,
                                                     torsion=torsion)

    def __str__(self):
        '''...'''
        s = super().__str__()
        s = s.replace(', ', ',')
        return s.replace('(', 'q^(')

    def _reduce_rep(self):
        '''reduces mod p the keys and values and deletes keys with 0 value
        or which are degenerate'''

        # reducing keys mod order
        if self.order != 'infinite':
            aux = list(self.items())
            self.clear()
            for x, v in aux:
                y = tuple(i % self.order for i in x)
                self[y] += v

        super()._reduce_rep()

    def phi(self):
        '''from Cyclic_DGModule to Barrat_Eccles'''

        if self.order == 'infinite':
            raise('phi is not define for elements of infinite order')

        r = self.order
        answer = {}
        for k, v in self.items():
            x = []
            for i in k:
                x.append(tuple(j % r + 1 for j in range(i, r + i)))
            answer[tuple(x)] = v

        return BarrattEccles_element(answer, torsion=self.torsion)

    def set_order(self, r):
        '''...'''
        setattr(self, 'order', r)
        self._reduce_rep()
        return self


class CubicalEilenbergZilber_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        '''...'''

        # check input data: dict with tuple of int keys
        if data:
            if not (isinstance(data, dict)
                    and all(isinstance(x, tuple) for x in data.keys())
                    and all(isinstance(i, str) for i in
                            chain.from_iterable(data.keys()))
                    ):
                raise TypeError(
                    'data type must be dict with tuple of str keys')

        # initialize element
        super(type(self), self).__init__(data=data, torsion=torsion)

    @ property
    def arity(self):
        if not self:
            return None

        arities = set(len(x) for x in self.keys())

        if len(arities) > 1:
            return arities

        return arities.pop()

    def coproduct(self):
        def _coproduct_mod2(x, answer=set()):
            '''mod 2 version'''
            if len(x) == 0:
                return {('', '')}

            if x[0] == 'e':
                answer = {('0' + pair[0], 'e' + pair[1])
                          for pair in _coproduct_mod2(x[1:], answer)} \
                    ^ {('e' + pair[0], '1' + pair[1])
                       for pair in _coproduct_mod2(x[1:], answer)}
            else:
                answer = {(x[0] + pair[0], x[0] + pair[1])
                          for pair in _coproduct_mod2(x[1:], answer)}
            return answer

        def _sign(pair):
            num_transps = 0
            for i, x in enumerate(pair[0]):
                if x == 'e':
                    num_transps += len([y for y in pair[1][:i + 1] if y == 'e'])
            return (-1)**num_transps

        answer = self.zero()
        if self.arity == 1:
            for k, v in self.items():
                for pair in _coproduct_mod2(k[0]):
                    answer += type(self)(
                        {pair: _sign(pair) * v}, torsion=self.torsion)
            return answer
        else:
            for k, v in self.items():
                x = type(self)(
                    {(k[0],): v}, torsion=self.torsion).coproduct()
                answer += type(self)(
                    {m + k[1:]: w for m, w in x.items()}, torsion=self.torsion)
            return answer

    def product(self, other=None):
        '''...'''

        def _product(self, other):
            '''...'''
            assert self.torsion == other.torsion, 'defined for equal torsion'
            assert self.arity == other.arity == 1, 'defined for arity 1'

            answer = other.zero()
            if not self or not other:
                return answer

            _ast = {'01': 'e', '10': 'e'}
            for k1, v1 in self.items():
                for k2, v2 in other.items():
                    x, y = k1[0], k2[0]
                    expected_degree = x.count('e') + y.count('e') + 1
                    for i, pair in enumerate(zip(x, y)):
                        try:
                            new_summand = y[:i] + _ast[pair[0] + pair[1]] + x[i + 1:]
                            if new_summand.count('e') == expected_degree:
                                answer += self.create({
                                    (new_summand,): v1 * v2 * (-1)**x.count('e')})
                        except KeyError:
                            pass
            return answer

        if other:
            return _product(self, other)

        answer = self.zero()
        for k, v in self.items():
            elements = [self.create({(a,): 1}) for a in k]
            answer += v * reduce(lambda x, y: x.product(y), elements)
        return answer

    def __call__(self, other):
        '''...'''
        assert self.torsion == other.torsion, "defined for the same ring"
        assert isinstance(self.degree, int) and self.arity, "defined for homogeneous surjections"
        if isinstance(other, CubicalEilenbergZilber_element):
            assert other.arity == 1, "defined for chains on a single cube"
            answer = other.zero()
            iterated_diagonal = other
            for _ in range(self.degree + self.arity - 1):
                iterated_diagonal = iterated_diagonal.coproduct()
            # print(iterated_diagonal)
            for k1, v1 in self.items():
                for k2, v2 in iterated_diagonal.items():
                    # sign
                    odds = [i for i, x in enumerate(k2) if x.count('e') % 2]
                    sign = v1 * v2
                    for idx, i in enumerate(odds):
                        sign *= (-1)**len([j for j in odds[idx + 1:] if k1[i] > k1[j]])
                    # elements
                    elements = []
                    for s in range(1, max(k1) + 1):
                        element = other.create({tuple(k2[i] for i, s_i in enumerate(k1)
                                                      if s_i == s): 1})
                        elements.append(element.product())
                    if all(elements):
                        for multipair in product(*(element.items() for element in elements)):
                            new_key = tuple(pair[0][0] for pair in multipair)
                            coeff = sign * reduce(lambda i, j: i * j, (pair[1] for pair in multipair))
                            answer += answer.create({new_key: coeff})
            return answer


class CubicalEilenbergZilber_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        '''...'''

        # check input data: dict with tuple of int keys
        if data:
            if not (isinstance(data, dict)
                    and all(isinstance(x, tuple) for x in data.keys())
                    and all(isinstance(i, str) for i in
                            chain.from_iterable(data.keys()))
                    ):
                raise TypeError(
                    'data type must be dict with tuple of str keys')

        # initialize element
        super(type(self), self).__init__(data=data, torsion=torsion)

    @ property
    def arity(self):
        if not self:
            return None

        arities = set(len(x) for x in self.keys())

        if len(arities) > 1:
            return arities

        return arities.pop()

    def coproduct(self):
        def _coproduct_mod2(x, answer=set()):
            '''mod 2 version'''
            if len(x) == 0:
                return {('', '')}

            if x[0] == 'e':
                answer = {('0' + pair[0], 'e' + pair[1])
                          for pair in _coproduct_mod2(x[1:], answer)} \
                    ^ {('e' + pair[0], '1' + pair[1])
                       for pair in _coproduct_mod2(x[1:], answer)}
            else:
                answer = {(x[0] + pair[0], x[0] + pair[1])
                          for pair in _coproduct_mod2(x[1:], answer)}
            return answer

        def _sign(pair):
            num_transps = 0
            for i, x in enumerate(pair[0]):
                if x == 'e':
                    num_transps += len([y for y in pair[1][:i + 1] if y == 'e'])
            return (-1)**num_transps

        answer = self.zero()
        if self.arity == 1:
            for k, v in self.items():
                for pair in _coproduct_mod2(k[0]):
                    answer += type(self)(
                        {pair: _sign(pair) * v}, torsion=self.torsion)
            return answer
        else:
            for k, v in self.items():
                x = type(self)(
                    {(k[0],): v}, torsion=self.torsion).coproduct()
                answer += type(self)(
                    {m + k[1:]: w for m, w in x.items()}, torsion=self.torsion)
            return answer

    def product(self, other=None):
        '''...'''

        def _product(self, other):
            '''...'''
            assert self.torsion == other.torsion, 'defined for equal torsion'
            assert self.arity == other.arity == 1, 'defined for arity 1'

            answer = other.zero()
            if not self or not other:
                return answer

            _ast = {'01': 'e', '10': 'e'}
            for k1, v1 in self.items():
                for k2, v2 in other.items():
                    x, y = k1[0], k2[0]
                    expected_degree = x.count('e') + y.count('e') + 1
                    for i, pair in enumerate(zip(x, y)):
                        try:
                            new_summand = y[:i] + _ast[pair[0] + pair[1]] + x[i + 1:]
                            if new_summand.count('e') == expected_degree:
                                answer += self.create({
                                    (new_summand,): v1 * v2 * (-1)**x.count('e')})
                        except KeyError:
                            pass
            return answer

        if other:
            return _product(self, other)

        answer = self.zero()
        for k, v in self.items():
            elements = [self.create({(a,): 1}) for a in k]
            answer += v * reduce(lambda x, y: x.product(y), elements)
        return answer


class EilenbergZilber_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        '''...'''

        # checking input
        if data:
            # check input data: dict of tuple of tuple of int
            if not (all((isinstance(multiop, tuple)) for multiop in
                        data.keys())
                    and all((isinstance(op, tuple) for op in
                             chain.from_iterable(data.keys())))
                    and all((isinstance(i, int) for i in
                             chain.from_iterable(
                             chain.from_iterable(data.keys()))))
                    ):
                raise TypeError('keys must be tuple of tuple of int')

        # initialize object
        super(EilenbergZilber_element, self).__init__(data=data,
                                                      torsion=torsion)

    @ property
    def arity(self):
        arities = set(len(multiop) for multiop in self.keys())
        if len(arities) != 1:
            raise TypeError('keys must have same arity')
        return arities.pop()

    def __str__(self):
        '''...'''
        if not self:
            return '0'

        string = ''
        for multiop, coeff in self.items():
            if coeff != 1:
                string += str(coeff)
            string += '('
            for op in multiop:
                if not op:
                    d = 'id'
                else:
                    d = f'd_{"d_".join(str(i) for i in op)}'

                string += d + ')x('
            string = string[:-2] + ' + '
        return string[:-3]

    def __call__(self, other):
        '''...'''
        if isinstance(other, int):
            other = Module_element(
                {tuple(range(other + 1)): 1}, torsion=self.torsion)
        answer = Module_element(torsion=self.torsion)
        for k1, v1 in self.items():
            for k2, v2 in other.items():
                indices = range(len(k2))
                new_key = tuple(
                    tuple(k2[i] for i in indices if i not in op) for op in k1)
                answer += Module_element({new_key: v1 * v2},
                                         torsion=self.torsion)
        return answer

    def _reduce_rep(self):
        '''...'''

        # order face maps in increasing value
        self_data = dict(self)
        self.clear()
        for multiop, coeff in self_data.items():
            ordered_multiop = tuple()
            for op in multiop:
                ordered_op = EilenbergZilber_element._face_maps_sort(op)
                ordered_multiop += (ordered_op,)
            self[ordered_multiop] += coeff

        super()._reduce_rep()

    @ staticmethod
    def _face_maps_sort(face_maps):
        '''puts the face maps in canonical order d < ... < d using the
        simplicial identity d_i d_j = d_j d_{i+1} if i >= j'''

        face_maps = list(face_maps)
        for index in range(1, len(face_maps)):

            currentvalue = face_maps[index]
            position = index

            while (position > 0
                   and face_maps[position - 1] >= currentvalue):

                face_maps[position] = face_maps[position - 1] + 1
                position = position - 1

            face_maps[position] = currentvalue

        return tuple(face_maps)


class SteenrodProduct():
    '''...'''

    def __init__(self, r, i, torsion='free'):
        # checking input
        if not isinstance(r, int) and isinstance(i, int):
            raise TypeError('initialize with two int')

        # setting attributes read arity and degree
        self.r, self.i, self.torsion = r, i, torsion

    def as_CyclicDGMolule_element(self):
        '''...'''
        # generator of W
        e = CyclicModule_element({0: 1}, torsion=self.torsion, order=self.r)
        if self.i < 0:  # set to e = 0
            e = CyclicModule_element(torsion=self.torsion, order=self.r)

        return e.psi(self.i)

    def as_BarrattEccles_element(self):
        '''...'''
        return self.as_CyclicDGMolule_element().phi()

    def as_Surjection_element(self):
        '''...'''
        return self.as_BarrattEccles_element().table_reduction()


class SteenrodOperation(object):
    '''Models a chain level representative of P_s or bP_s over the prime p
    acting on an element of degree n'''

    def __init__(self, p, s, n, bockstein=False, convention='chain'):

        # check input
        if not (isinstance(p, int)
                and isinstance(s, int) and isinstance(n, int)):
            raise TypeError('initialize with three int p,s,n')
        if p == 2 and bockstein:
            raise TypeError('bP only defined for odd primes')
        if not isinstance(bockstein, bool):
            raise TypeError('bockstein must be a boolean')
        if convention != 'chain' and convention != 'cochain':
            raise TypeError("convention must be either 'chain' or 'cochain'")

        # setting attributes
        self.p, self.s, self.n = p, s, n
        self.b = int(bockstein)

        if convention == 'chain':
            self.c = 1
        elif convention == 'cochain':
            self.c = -1

        if p == 2:
            setattr(self, 'coeff', 1)
            # chain: s-n & cochain: -n-s
            setattr(self, 'd', (self.c) * (s - n))

        elif p > 2:
            # Serre convention: v(2j)=(-1)^j & v(2j+1)=v(2j)*m! w/ m=(p-1)/2
            coeff = (-1)**(floor(n / 2) + s)
            if n / 2 - floor(n / 2):
                coeff *= factorial((p - 1) / 2)
            setattr(self, 'coeff', int(coeff))
            # degree of e: chain (2s-n)(p-1)-b & cochain (n+2s)(p-1)-b
            setattr(self, 'd', (self.c) * (2 * s - n) * (p - 1) - (self.b))

    def __str__(self):
        string = f'(P^{self.s})_{{{self.n}}}'
        if self.b:
            string = string.replace('(', '(b')
        return string

    def as_CyclicDGMolule_element(self):
        '''...'''
        # generator of W
        e = CyclicModule_element({0: 1}, torsion=self.p, order=self.p)
        if self.d < 0:  # set to e = 0
            e = CyclicModule_element(torsion=self.p, order=self.p)

        return (self.coeff) * e.psi(self.d)

    def as_BarrattEccles_element(self):
        '''...'''
        return self.as_CyclicDGMolule_element().phi()

    def as_Surjection_element(self):
        '''...'''
        return self.as_BarrattEccles_element().table_reduction()

    def as_EilenbergZilber_element(self):
        '''...'''
        if self.c != -1:
            raise TypeError('EilenbergZilber uses cochain convention')
        if self.p == 2:
            return self.as_Surjection_element().interval_cut(
                self.n + self.s)
        if self.p > 2:
            return self.as_Surjection_element().interval_cut(
                self.n + (2 * self.s) * (self.p - 1) + self.b)
