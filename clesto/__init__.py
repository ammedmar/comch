from collections import Counter
from itertools import chain, combinations, permutations, tee, product
from functools import reduce
from math import floor, factorial
from operator import attrgetter

# Check compositions for finite torsion

# To-do
# Composition of EZ_elements
# Improve str functions
# Cell for composition of Permutations
# Cell for composition of BE_element
# Cell for cut_interval
# Revisit definition of EZ_elmt using n as a parameter


def partitions(n, k, smallest_value=1, largest_value=None, ordered=False):
    '''n is the integer to partition and k is the length of partitions.
    It returns all k tuples of integers greater or equal to smallest_value
    and less than or equal to largest_value that add up to n.
    If ordered == True it returns all tuples if False it returns those
    in non-decreassing order '''
    if largest_value is None:
        largest_value = n

    def unordered_partitions(n, k, r=smallest_value, m=largest_value):
        if k == 1:
            if r <= n <= m:
                yield (n,)
            return
        for i in range(r, m + 1):
            for result in unordered_partitions(n - i, k - 1, i, m):
                yield (i,) + result

    if ordered:
        return chain.from_iterable(set(permutations(p)) for p
                                   in unordered_partitions(n, k))
    if not ordered:
        return unordered_partitions(n, k)


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class Module_element(Counter):
    """
    Counter with arithmetic improvements to handle (modular) integer values.

    Class constructed to model free module elements over the ring of
    (modular) integers.

    Attributes
    ----------
    default_torsion : non-negative int or string 'free'
        An int m sets R = Z/mZ whereas 'free' sets R = Z

    """

    default_torsion = 'free'

    def __init__(self, data=None, torsion=None):
        # print('initializing as Module_element')

        # check input data: dict with int values
        if data:
            if not (isinstance(data, dict) and
                    all((type(v) is int for v in data.values()))
                    ):
                raise TypeError('input must be dict with int values')

        # checking input torsion: positive int or 'free'
        if torsion is not None:
            if not (isinstance(torsion, int) and torsion > 0 or
                    torsion == 'free'
                    ):
                raise TypeError("torsion must be a positive int or 'free'")

        # setting torsion
        m = torsion if torsion else type(self).default_torsion
        setattr(self, 'torsion', m)

        # initialize element
        super(Module_element, self).__init__(data)

        self._reduce_rep()

    def __str__(self):
        if not self:
            return '0'
        else:
            answer = ''
            for key, value in self.items():
                if value < -1:
                    answer += f'- {abs(value)}{key} '
                elif value == -1:
                    answer += f'- {key} '
                elif value == 1:
                    answer += f'+ {key} '
                elif value > 1:
                    answer += f'+ {value}{key} '
            if answer[0] == '+':
                answer = answer[2:]

            return answer[:-1]

    def __add__(self, other):
        '''The sum of two free module elements.

        >>> Module_element({'a':1, 'b':2}) + Module_element({'a':1})
        Module_element({'a':2, 'b':2})

        '''
        self.compare_attributes(other)
        answer = type(self)(self).copy_attrs_from(self)
        answer.update(other)
        answer._reduce_rep()
        return answer

    def __sub__(self, other):
        '''The substraction of two free module elements.

        >>> Module_element({'a':1, 'b':2}) - Module_element({'a':1})
        Module_element({'b':2})

        '''
        self.compare_attributes(other)
        answer = type(self)(self).copy_attrs_from(self)
        answer.subtract(other)
        answer._reduce_rep()
        return answer

    def __rmul__(self, c):
        '''The scaling by c of a free module element.

        >>> 3*Module_element({'a':1, 'b':2})
        Module_element({'a':3, 'b':6})

        '''
        if not isinstance(c, int):
            raise TypeError(f'can not act by non-int of type {type(c)}')

        scaled = {k: c * v for k, v in self.items()}
        return type(self)(scaled).copy_attrs_from(self)

    def __neg__(self):
        '''The additive inverse of a free module element.

        >>> -Module_element({'a':1, 'b':2})
        Module_element({'a':-1, 'b':-22})

        '''
        return self.__rmul__(-1)

    def __iadd__(self, other):
        '''The in place addition of two free module elements.

        >>> x = Module_element({'a':1, 'b':2})
        >>> x += Module_element({'a':3, 'b':6})
        Module_element({'a':4, 'b':8})

        '''
        self.compare_attributes(other)
        self.update(other)
        self._reduce_rep()
        return self

    def __isub__(self, other):
        '''The in place addition of two free module elements.

        >>> x = Module_element({'a':1, 'b':2})
        >>> x -= Module_element({'a':3, 'b':6})
        Module_element({'a':-2, 'b':-4})

        '''
        self.compare_attributes(other)
        self.subtract(other)
        self._reduce_rep()
        return self

    def _reduce_rep(self):
        '''The preferred representative of the free module element.

        It reduces all values mod n if torsion is n and removes
        key:value pairs with value = 0.

        >>> Module_element({'a':1, 'b':2, 'c':0})
        Module_element({'a':1, 'b':2})

        '''
        # print('reducing as Module_element')
        # reducing coefficients mod torsion
        if self.torsion != 'free':
            for key, value in self.items():
                self[key] = value % self.torsion

        # removing key:value pairs with value = 0
        zeros = [k for k, v in self.items() if not v]
        for key in zeros:
            del self[key]

    def set_torsion(self, m):
        '''...'''
        setattr(self, 'torsion', m)
        self._reduce_rep()
        return self

    def compare_attributes(self, other):
        '''...'''
        if self.__dict__ != other.__dict__:
            raise AttributeError('not the same attributes')

    def copy_attrs_from(self, other):
        '''...'''
        for attr, value in other.__dict__.items():
            setattr(self, attr, value)
        self._reduce_rep()
        return(self)


class CyclicModule_element(Module_element):
    '''Modeling elements in Z/mZ[C_n]'''

    default_order = 'infinite'

    def __init__(self, data=None, torsion=None, order=None):
        # print("initializing as CyclicModule_element")

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
        # print('reducing as CyclicModule_element')

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
        '''Uses the direct formula not the one involving the table_reduction morphism'''

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
        coeff = reduce(lambda x, y: x+y, rhos).set_torsion(m) # sum of all rhos

        return coeff * Surjection_element({key: 1 for key in keys}, torsion=m)

class SymmetricModule_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        # print("initializing as SymmetricModule_element")

        # check input data: dict with tuple of int keys
        if data:
            if not (isinstance(data, dict) and
                    all(isinstance(perm, tuple) for perm in data.keys()) and
                    all(isinstance(i, int) for i in
                        chain.from_iterable(data.keys()))
                    ):
                raise TypeError('data type must be dict with tuple of int keys')

            if any((set(k) != set(range(1, len(k) + 1)) for k in data.keys())):
                raise TypeError('keys must be permutations of (1,2,...,r)')

        # initialize element
        super(SymmetricModule_element, self).__init__(data=data,
                                                      torsion=torsion)

    @property
    def arity(self):
        '''...'''
        if not self:
            return None
        
        arities = set(max(k) for k in self.keys())
        if len(arities) > 1:
            return arities

        return arities.pop()

    def __str__(self):
        '''...'''
        if not self:
            return '0'
        else:
            s = super().__str__()
            return s.replace(', ', ',')

    def __mul__(self, other):
        '''...'''
        if isinstance(other, SymmetricModule_element):
            self.compare_attributes(other)
            answer = type(other)().copy_attrs_from(self)
            for k1, v1 in self.items():
                for k2, v2 in other.items():
                    answer[tuple(k1[i - 1] for i in k2)] += v1 * v2
            answer._reduce_rep()
            return answer

        if isinstance(other, BarrattEccles_element):
            self.compare_attributes(other)
            answer = type(other)().copy_attrs_from(other)
            for k1, v1 in self.items():
                for k2, v2 in other.items():
                    new_key = tuple(tuple(k1[i - 1] for i in x) for x in k2)
                    answer[new_key] = v1 * v2
            answer._reduce_rep()
            return answer

        if isinstance(other, Surjection_element):
            self.compare_attributes(other)
            answer = type(other)().copy_attrs_from(other)
            for k1, v1 in self.items():
                for k2, v2 in other.items():
                    new_key = tuple(k1[i - 1] for i in k2)
                    answer[new_key] = v1 * v2
            answer._reduce_rep()
            return answer

    def __rmul__(self, other):
        '''scaled by integers on the left and acts on the right'''

        if isinstance(other, int):
            return super().__rmul__(other)

        if isinstance(other, BarrattEccles_element):
            self.compare_attributes(other)
            answer = type(other)().copy_attrs_from(other)
            for k1, v1 in self.items():
                for k2, v2 in other.items():
                    new_key = tuple(tuple(x[i - 1] for i in k1) for x in k2)
                    answer[new_key] = v1 * v2
            answer._reduce_rep()
            return answer

    def __call__(self, other):
        '''...'''
        return self.__mul__(other)

    def compose(self, *others):
        '''...'''
        # set compose with 0 equal 0
        if self == SymmetricModule_element():
            return SymmetricModule_element()

        # partial composition
        if len(others) == 2 and isinstance(others[1], int):
            # unpack and check data
            other, k = others
            if self.torsion != other.torsion:
                raise TypeError('elements must have equal torsion')

            # initialize answer
            comp = SymmetricModule_element(torsion=self.torsion)

            # populate answer using linearity
            for perm1, coeff1 in self.items():
                for perm2, coeff2 in other.items():
                    s = len(perm2) - 1
                    to_insert = tuple(i + k - 1 for i in perm2)
                    at = perm1.index(k)
                    shift = tuple(map(lambda i: i + s if i > k else i, perm1))
                    inserted = shift[:at] + to_insert + shift[at + 1:]
                    new_coeff = coeff1 * coeff2
                    to_add = SymmetricModule_element({inserted: new_coeff},
                                                     torsion=self.torsion)
                    comp += to_add
            return comp

        # total composition
        else:
            if len(others) != self.arity:
                raise TypeError('argument number must equal the arity of' +
                                f'self, but {len(others)} != {self.arity}')
            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer.compose(other, idx + 1)
            return answer

    @classmethod
    def rho(self, arity, exponent=1, torsion='free'):
        key = tuple(i % arity + 1 for i in range(exponent, exponent + arity))
        return SymmetricModule_element({key: 1}, torsion=torsion)

class DGModule_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        # print('initializing as DGModule_element')

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
        # print('reducing as DGModule_element')

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
        # print("initializing as CyclicDGModule_element")

        # check input data: dict with tuple of int keys
        if data:
            if not (isinstance(data, dict) and
                    all((isinstance(i, int) for i in
                         chain.from_iterable(data.keys())))
                    ):
                raise ValueError('data type must be dict' +
                                 'with tuple of int keys')
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
        # print('reducing as CyclicDGModule_element')

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

class BarrattEccles_element(DGModule_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        # print("initializing as BarrattEccles_element")

        # check input data: dict with tuple of tuple of int keys
        if data:
            if not (isinstance(data, dict) and
                    all(isinstance(x, tuple) for x in data.keys()) and
                    all(isinstance(perm, tuple) for perm in
                        chain.from_iterable(data.keys())) and
                    all(isinstance(i, int) for i in
                        chain.from_iterable(chain.from_iterable(data.keys())))
                    ):
                raise TypeError('data type must be dict ' +
                                'with tuple of tuple of int keys')

            if any((set(perm) != set(range(1, len(perm) + 1)) for perm in
                    chain.from_iterable(data.keys()))):
                raise TypeError('keys must tuples of ' +
                                'permutations of (1,2,...,r)')

        # initializing element
        super(BarrattEccles_element, self).__init__(data=data,
                                                    torsion=torsion)
    @property
    def arity(self):
        '''...'''
        if not self:
            return None

        arities = set()
        for k in self.keys():
            arities_in_k = set(max(perm) for perm in k)
            arities |= arities_in_k  # in place union

        if len(arities) > 1:
            return arities

        return arities.pop()

    def _paths(p, q):
        '''returns as a list all increasing paths from (0,0) to (p,q)'''

        if (p, q) == (0, 0):
            return [((0, 0),)]

        answer = list()
        if p > 0:
            west = BarrattEccles_element._paths(p - 1, q)
            for path in west:
                answer.append(path + ((p, q),))

        if q > 0:
            south = BarrattEccles_element._paths(p, q - 1)
            for path in south:
                answer.append(path + ((p, q),))

        return answer

    def _sgn_of_path(path):
        '''...'''
        segments = range(1, len(path))
        horizontal_segments = []
        vertical_segments = []
        for i in segments:
            vertex1 = path[i - 1]
            vertex2 = path[i]
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
                sgn = diff // abs(diff)
        return sgn

    def compose(self, *others):
        '''...'''

        # partial composition
        if len(others) == 2 and isinstance(others[1], int):
            # unpaking and checking input
            other, k = others
            if self.torsion != other.torsion:
                raise TypeError('not the same torsion')

            # initialize answer
            answer = BarrattEccles_element(torsion=self.torsion)

            # populate answer using linearity
            for perm_vect1, coeff1 in self.items():
                for perm_vect2, coeff2 in other.items():
                    comp = BarrattEccles_element().copy_attrs_from(answer)
                    p, q = len(perm_vect1) - 1, len(perm_vect2) - 1
                    # summands parametrized by paths from (0,0) to (p,q)
                    for path in BarrattEccles_element._paths(p, q):
                        new_perm_vect = ()
                        for i, j in path:
                            perm1 = SymmetricModule_element({perm_vect1[i]: 1}, torsion=self.torsion)
                            perm2 = SymmetricModule_element({perm_vect2[j]: 1}, torsion=self.torsion)
                            partial_comp = perm1.compose(perm2, k)
                            new_perm_vect += (tuple(partial_comp.keys())[0],)
                        sgn = BarrattEccles_element._sgn_of_path(path)
                        comp += BarrattEccles_element({new_perm_vect: sgn}, torsion=self.torsion)
                    answer += coeff1 * coeff2 * comp
            return answer

        # total composition
        else:
            if not len(others) == self.arity:
                raise TypeError('the number of arguments must be equal to ' +
                                'the arity of self')
            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer.compose(other, idx + 1)
            return answer

    def table_reduction(self):
        '''given a set of basis element in the Barratt_Eccles operad, it returns
        the set of surjections in its image via the table reduction morphism'''

        answer = Surjection_element(torsion=self.torsion)

        for bar_ecc_element, value in self.items():
            d, a = len(bar_ecc_element) - 1, max(bar_ecc_element[0])
            for pi in partitions(d + a, d + 1, ordered=True):
                surjection, removed = [], []
                degenerate = False
                for idx, i in enumerate(pi):
                    filtered = [i for i in bar_ecc_element[idx]
                                if i not in removed]
                    if idx > 0 and surjection[-1] == filtered[0]:
                        degenerate = True
                        break
                    if i > 1:
                        removed += filtered[: i - 1]
                    surjection += filtered[: i]

                if not degenerate:
                    answer += Surjection_element({tuple(surjection): value},
                                                 torsion=self.torsion)
        answer._reduce_rep()
        return answer

class Surjection_element(DGModule_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        '''...'''
        # print('initializing as Surjection_element')
        # check input data: dict with tuple of int keys
        if data:
            if not (isinstance(data, dict) and
                    all(isinstance(surj, tuple) for surj in data.keys()) and
                    all(isinstance(i, int) for i in
                        chain.from_iterable(data.keys()))
                    ):
                raise TypeError('data type must be dict with tuple of int keys')

        # initialize element
        super(Surjection_element, self).__init__(data=data, torsion=torsion)

    def _reduce_rep(self):
        '''...'''
        # print('reducing as Surjection_element')
        zeros = list()
        for surj in self.keys():
            if set(surj) != set(range(1, max(surj) + 1)):
                zeros.append(surj)
        for k in zeros:
            del self[k]

        super()._reduce_rep()

    @property
    def arity(self):
        if not self:
            return None

        arities = set(max(surj) for surj in self.keys())

        if len(arities) > 1:
            return arities
        
        return arities.pop()
    
    @property
    def complexity(self):
        '''returns the complexity of an element in the Surjection operad'''
        complexity = {1}
        for surjection in self.keys():
            for i, j in combinations(range(1, max(surjection)+1), 2):
                r = tuple(k for k in surjection if k == i or k == j)
                if all([i != j for i, j in pairwise(r)]) and len(r) > 1:
                    complexity.add(len(r)-1)
        return max(complexity)

    def table_arrangement(surj, only_dict=False):
        '''Returns the table arrangement of a surjection, as a tuple.
        If only_dict=True, it returns a dictionary dict such that
        dict[index] is the table arrangement of surj where
        surj[index] lies.
        '''
        finals = set()
        for k in range(1, max(surj) + 1):
            for index in reversed(range(len(surj))):
                if surj[index] == k:
                    finals.add(index)
                    break

        if only_dict:
            result = dict()
            row = 0
            for index, el in enumerate(surj):
                result[index] = row
                if index not in finals or index == len(surj) - 1:
                    row += 1
            return result

        result = tuple()
        new_row = tuple()
        for index, el in enumerate(surj):
            new_row += (el,)
            if index not in finals or index == len(surj) - 1:
                result += (new_row,)
                new_row = tuple()
        return result

    def last_index(value, sequence):
        '''returns the index of the last occurence of a value in a sequence'''

        for index, element in reversed(list(enumerate(sequence))):
            if value == element:
                return index
        raise ValueError(f"{value} is not an element of {sequence}")

    def degrees(surj, pieces):
        '''returns a tuple with the degree of the pieces in surj'''

        row = Surjection_element.table_arrangement(surj, only_dict=True)
        result = tuple()
        remaining_surj = surj
        for piece in reversed(pieces):
            end = piece[-1]
            end_index = Surjection_element.last_index(end, remaining_surj)
            start_index = end_index - len(piece) + 1
            degree = row[end_index] - row[start_index]
            result = (degree,) + result
            remaining_surj = remaining_surj[0:start_index + 1]
        return result

    def _pcompose(u, v, k):
        '''partial composition of Surjection_elements at place k'''

        result = Surjection_element(torsion=u.torsion)
        for surj1, coeff1 in u.items():
            for surj2, coeff2 in v.items():
                # shift surj1 and surj2
                surj1_shifted = tuple(map(lambda i: i + v.arity - 1
                                          if i > k else i, surj1))
                surj2_shifted = tuple(el + k - 1 for el in surj2)

                # find the occurences of k in surj1
                occurences = tuple(index for index, el in enumerate(surj1)
                                   if el == k)

                # cut surj1 at the positions contained in occurences
                indices_to_cut_surj1 = (0,) + occurences + (len(surj1) - 1,)

                surj1_cut = tuple()
                for i in range(len(indices_to_cut_surj1) - 1):
                    left = indices_to_cut_surj1[i]
                    right = indices_to_cut_surj1[i + 1]
                    surj1_cut += (surj1_shifted[left:right + 1],)

                # cut surj2 into len(occurences) - 1 pieces in
                # all possible ways
                num_pieces = len(occurences) - 1
                cuts_for_surj2 = combinations(range(len(surj2)), num_pieces)
                for indices_to_cut_surj2 in cuts_for_surj2:
                    surj2_cut = tuple()
                    inserted = surj1_shifted
                    indices_to_cut_surj2 = (0,) + indices_to_cut_surj2
                    indices_to_cut_surj2 += (len(surj2) - 1,)
                    for i in reversed(range(len(indices_to_cut_surj2) - 1)):
                        left = indices_to_cut_surj2[i]
                        right = indices_to_cut_surj2[i + 1]
                        to_insert = surj2_shifted[left:right + 1]
                        surj2_cut = (to_insert,) + surj2_cut
                        index = occurences[i]
                        left_inserted = inserted[0:index]
                        right_inserted = inserted[index + 1:len(inserted)]
                        inserted = left_inserted + to_insert + right_inserted

                    # compute degrees of pieces of the cut of surj1
                    deg1 = Surjection_element.degrees(surj1_shifted, surj1_cut)
                    deg2 = Surjection_element.degrees(surj2_shifted, surj2_cut)

                    # compute the sign of the permutation
                    sign_exp = 0
                    for index2 in range(len(surj2_cut)):
                        deg_piece2 = deg2[index2]
                        left = occurences[index2]
                        right = len(surj1_cut)
                        deg_pieces1 = sum(deg1[index1]
                                          for index1 in range(left, right))
                        sign_exp += (deg_piece2 * deg_pieces1) % 2
                    sign = (-1)**sign_exp

                    # add the new element to the result
                    new_coeff = coeff1 * coeff2 * sign
                    result += Surjection_element({inserted: new_coeff}, torsion=v.torsion)
        return result

    def compose(self, *others):
        '''general composition of Surjection_elements'''

        # partial composition
        if len(others) == 2 and isinstance(others[1], int):
            # unpaking and checking input
            other, k = others
            if self.torsion != other.torsion:
                raise TypeError('not the same torsion')
            return self._pcompose(other, k)

        # total composition
        else:
            if not len(others) == self.arity:
                raise TypeError('the number of arguments must be equal to ' +
                                'the arity of self')
            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer._pcompose(other, idx + 1)
            return answer

    def interval_cut(self, n):
        '''...'''
        def all_cuts(k, n):
            '''(k,n) -> 0 = n0 <= n1 <= ... <= nk <= n'''
            if k == 0:
                yield (0,)
            if k > 0:
                for cut in all_cuts(k - 1, n):
                    for i in range(cut[-1], n + 1):
                        yield cut + (i,)

        def constraints(surj):
            '''...'''
            for i in range(1, max(surj) + 1):
                preimage = (idx for idx, s in enumerate(surj) if s == i)
                for s, t in pairwise(preimage):
                    yield (s + 1, t)

        def good_cut(cut, const):
            '''...'''
            return all(cut[i] != cut[ipp] for i, ipp in const)

        def cut2multioperator(cut, n):
            '''...'''
            multisimplex = {i: tuple() for i in range(1, max(surj) + 1)}
            for i, pair in enumerate(pairwise(cut)):
                multisimplex[surj[i]] += tuple(range(pair[0], pair[1] + 1))

            std = set(range(n + 1))  # {0,1,...,n}
            return tuple(tuple(std.difference(set(spx)))
                         for spx in multisimplex.values())

        def cut_sign(cut, surj):

            class LeveledInterval:
                def __init__(self, start, end, level):
                    self.start = start
                    self.end = end
                    self.level = level
                    self.is_inner = True

                def length(self):
                    length = interval.end - interval.start
                    if self.is_inner:
                        length += 1
                    return length

            # transform cut to tuple of intervals
            intervals = tuple(LeveledInterval(cut[i], cut[i + 1], surj[i])
                              for i in range(len(surj)))

            # classify intervals: internal & final
            for level in range(1, max(surj) + 1):
                for interval in reversed(intervals):
                    if level == interval.level:
                        interval.is_inner = False
                        break

            # position sign = sum of n_i over internal intervals (n_{i-1}, n_i)
            position_sign_exp = sum(interval.end for interval in intervals
                                    if interval.is_inner)

            # permutation sign associated to ordering the surjection
            ordered = sorted(intervals, key=attrgetter('level'))
            perm_sign_exp = sum(ordered[i].length() * ordered[j].length()
                                for i in range(len(ordered))
                                for j in range(i, len(ordered))
                                if intervals.index(ordered[i]) >
                                intervals.index(ordered[j]))

            return (-1)**((position_sign_exp + perm_sign_exp) % 2)

        answer = EilenbergZilber_element().copy_attrs_from(self)
        for surj, coeff in self.items():
            const = set(constraints(surj))
            k = len(surj) - 1
            good_cuts = (cut + (n,) for cut in all_cuts(k, n)
                         if good_cut(cut + (n,), const))
            for cut in good_cuts:
                sign = cut_sign(cut, surj)
                multiop = cut2multioperator(cut, n)
                answer += EilenbergZilber_element(
                    {multiop: sign * coeff}).copy_attrs_from(answer)
        return answer


class EilenbergZilber_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        '''...'''
        # print('initializing as EilenbergZilber_element')
        # check input
        if data:
            # check input data: dict of tuple of tuple of int
            if not (all((isinstance(multiop, tuple)) for multiop in
                        data.keys()) and
                    all((isinstance(op, tuple) for op in
                         chain.from_iterable(data.keys()))) and
                    all((isinstance(i, int) for i in
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
            other = Module_element({tuple(range(other+1)):1}, torsion=self.torsion)
        answer = Module_element(torsion=self.torsion)
        for k1, v1 in self.items():
            for k2, v2 in other.items():
                indices = range(len(k2))
                new_key = tuple(tuple(k2[i] for i in indices if not i in op) for op in k1)
                answer += Module_element({new_key: v1 * v2}, torsion=self.torsion)
        return answer

    def _reduce_rep(self):
        '''...'''
        # print('reducing as EilenbergZilber_element')

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

    @staticmethod
    def _face_maps_sort(face_maps):
        '''puts the face maps in canonical order d < ... < d using the
        simplicial identity d_i d_j = d_j d_{i+1} if i >= j'''

        face_maps = list(face_maps)
        for index in range(1, len(face_maps)):

            currentvalue = face_maps[index]
            position = index

            while (position > 0 and
                   face_maps[position - 1] >= currentvalue):

                face_maps[position] = face_maps[position - 1] + 1
                position = position - 1

            face_maps[position] = currentvalue

        return tuple(face_maps)


class SteenrodOperation(object):
    '''Models a chain level representative of P^s or bP^s over the prime p
    acting on an element of degree n'''

    def __init__(self, p, s, n, bockstein=False, convention='chain'):

        # check input
        if not (isinstance(p, int) and
                isinstance(s, int) and isinstance(n, int)):
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
