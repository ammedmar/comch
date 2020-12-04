from .module import Module_element, TorsionError
from itertools import chain, product


class ArityError(Exception):
    """Exception raised for unequal arity.

    Attributes:
        message : explanation of the error
    """

    def __init__(self, message='attribute arity must be the same'):
        self.message = message
        super(ArityError, self).__init__(message)


class SymmetricGroup_element(tuple):
    '''Tuple with arithmetic improvements modeling permutations in
    a finite symmetric group.

    '''

    def __init__(self, perm):
        '''...'''
        if set(perm) != set(range(1, len(perm) + 1)):
            raise TypeError('must be a permutations of (1,2,...,r)')

    def __str__(self):
        s = super().__str__()
        return s.replace(', ', ',')

    @property
    def sign(self):
        '''Returns the sign of the permutation.

        >>> SymmetricGroup_element((5, 2, 4, 3, 1)).sign
        1

        '''

        def to_cycles(self, singletons=False):
            '''Transforms from bijection to collection of cycles.

            >>> SymmetricGroup_element((5, 2, 4, 3, 1)).to_cycles()
            [(1, 5), (3, 4)]

            '''
            p = list(self)
            cycles = []
            for i in range(len(p)):
                if p[i] is False:
                    continue
                cycleFirst = i + 1
                cycle = [cycleFirst]
                p[i], next = False, p[i]
                while next != cycleFirst:
                    cycle.append(next)
                    p[next - 1], next = False, p[next - 1]
                # add the cycle to the list of cycles
                if singletons or len(cycle) > 1:
                    cycles.append(tuple(cycle))

            return cycles

        cycles = to_cycles(self)
        return (-1)**sum(len(cycle) - 1 for cycle in cycles)

    @property
    def arity(self):
        '''Returns the arity of the permutation, i.e., the cardinality of
        its domain thought of as a bijection.

        >>> SymmetricGroup_element((5, 2, 4, 3, 1)).arity
        5

        '''
        return max(self)

    def inverse(self):
        '''Returns the inverse of the given symmetric group element.

        >>> pi = SymmetricGroup_element((2,3,1))
        >>> print(pi.inverse())
        (3,1,2)

        '''
        inverse = tuple(self.index(i + 1) + 1 for i in range(self.arity))
        return SymmetricGroup_element(inverse)

    def __mul__(self, other):
        '''product in the symmetric group: self * other ordered as bijections.

        >>> x = SymmetricGroup_element((2, 3, 1))
        >>> y = SymmetricGroup_element((1, 3, 2))
        >>> print(x * y)
        (2,1,3)

        '''
        if not isinstance(other, SymmetricGroup_element):
            self = SymmetricRing_element({self: 1}, torsion=other.torsion)
            return other.__rmul__(self)

        if self.arity != other.arity:
            raise TypeError('must have the same arity')

        return SymmetricGroup_element(tuple(self[i - 1] for i in other))

    def __pow__(self, times):
        '''iterated product in the symmetric group: self * ... * self.

        >>> x = SymmetricGroup_element((2, 3, 4, 5, 1))
        >>> print(x**5)
        (1,2,3,4,5)

        '''
        if times == 0:
            return SymmetricGroup_element(tuple(range(1, max(self) + 1)))
        answer = self
        for i in range(times - 1):
            answer *= self

        return answer

    def compose(self, *others):
        '''Operadic composition; read as composition of functions
        x o (y_1, ..., y_r) from right to left.

        >>> x = SymmetricGroup_element((1, 3, 2))
        >>> y = SymmetricGroup_element((2, 1))
        >>> print(x.compose(y, 1))
        (2,1,4,3)

        '''
        # partial composition
        if len(others) == 2 and isinstance(others[1], int):
            other, j = others
            s = len(other) - 1
            to_insert = tuple(i + j - 1 for i in other)
            at = self.index(j)
            shift = tuple(map(lambda i: i + s if i > j else i, self))
            answer = shift[:at] + to_insert + shift[at + 1:]
            return SymmetricGroup_element(answer)

        # total composition
        else:
            if len(others) != self.arity:
                raise TypeError('argument number must equal the arity of self')
            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer.compose(other, idx + 1)
            return answer


class SymmetricRing():
    '''Class to produce special elements in the group ring of finite
    symmetric groups

    '''

    @staticmethod
    def identity_element(arity, torsion=None):
        identity = tuple(range(1, arity + 1))
        return SymmetricRing_element({identity: 1}, torsion=torsion)

    @staticmethod
    def rotation_element(arity, torsion=None):
        rho = tuple(range(2, arity + 1)) + (1,)
        return SymmetricRing_element({rho: 1}, torsion=torsion)

    @staticmethod
    def transposition_element(arity, torsion=None):
        rho = tuple(range(2, arity + 1)) + (1,)
        identity = tuple(range(1, arity + 1))
        return SymmetricRing_element({rho: 1, identity: -1},
                                     torsion=torsion)

    @staticmethod
    def norm_element(arity, torsion=None):
        rho = SymmetricRing.rotation_element(arity, torsion=torsion)
        answer = SymmetricRing_element(torsion=torsion)
        for i in range(arity):
            answer += rho**i
        return answer


class SymmetricRing_element(Module_element):
    '''Module_element modeling elements in the integral group ring of finite
    symmetric groups.

    '''

    def __init__(self, data=None, torsion=None):

        if data:
            if not (isinstance(data, dict)
                    and all(isinstance(perm, tuple) for perm in data.keys())
                    and all(isinstance(i, int) for i in
                            chain.from_iterable(data.keys()))
                    ):
                raise TypeError(
                    'data type must be dict with tuple of int keys')

            data = {SymmetricGroup_element(k): v for k, v in data.items()}

        super(SymmetricRing_element, self).__init__(data=data,
                                                    torsion=torsion)

    def __str__(self):
        s = super().__str__()
        return s.replace(', ', ',')

    @property
    def arity(self):
        ''' Returns the cardinality of the domain of the involved permutations

        >>> SymmetricRing_element({(5, 2, 4, 3, 1): 1}).arity
        5
        >>> SymmetricRing_element({(2, 3, 1): 1, (1, 2): 1}).arity

        '''
        if not self:
            return None

        arities = set(max(k) for k in self.keys())
        if len(arities) > 1:
            return None

        return arities.pop()

    def __mul__(self, other):
        ''' left multiplication by symmetric ring element

        >>> p = SymmetricRing_element({(4, 3, 2, 1): 1, (1, 2, 3, 4): 2})
        >>> print(3 * p)
        3(4,3,2,1) + 6(1,2,3,4)
        >>> q = SymmetricRing_element({(4, 1, 2, 3): 1})
        >>> print(p * q)
        (1,4,3,2) + 2(4,1,2,3)

        '''

        if isinstance(other, int):
            return super().__rmul__(other)

        if not isinstance(other, SymmetricRing_element):
            return other.__rmul__(self)

        if self.torsion != other.torsion:
            raise TorsionError

        answer = self.zero()
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            answer[tuple(k1[i - 1] for i in k2)] += v1 * v2

        answer._reduce_rep()
        return answer

    def __pow__(self, times):
        '''...'''

        if times == 0:
            return SymmetricRing.identity_element(self.arity, self.torsion)
        answer = self.zero()
        for k, v in self.items():
            answer += self.create({k**times: v})
        return answer

    def compose(self, *others):
        '''Operadic composition, read as composition of functions
        x o (y_1, ..., y_r) from right to left.

        >>> x = SymmetricRing_element({(2, 3, 1): 1, (1, 2, 3): -1})
        >>> y = SymmetricRing_element({(2, 1): 1, (1, 2): 1})
        >>> print(x.compose(y, 2))
        (3,2,4,1) + (2,3,4,1) - (1,3,2,4) - (1,2,3,4)

        '''

        if not self:
            return self.zero()

        # partial composition
        if len(others) == 2 and isinstance(others[1], int):

            other, j = others
            if self.torsion != other.torsion:
                raise TypeError('elements must have equal torsion')

            answer = self.zero()
            for (k1, v1), (k2, v2) in product(self.items(), other.items()):
                new_k = k1.compose(k2, j)
                new_v = v1 * v2
                to_add = self.create({new_k: new_v})
                answer += to_add
            return answer

        # total composition
        else:
            if len(others) != self.arity:
                raise TypeError('argument number must equal arity of self')
            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer.compose(other, idx + 1)
            return answer
