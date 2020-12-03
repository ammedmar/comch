from ..basics import Module_element, TorsionError
from ..basics import SymmetricGroup_element, ArityError
from ..basics import SymmetricModule_element, SymmetricModule

from ..eilenberg_zilber import Simplex, EilenbergZilber_element, EilenbergZilber
from ..eilenberg_zilber import CubicalEilenbergZilber_element, CubicalEilenbergZilber
from ..utils import pairwise

from itertools import chain, combinations, product
from operator import itemgetter
from functools import reduce
from math import floor, factorial


class Surjection_element(Module_element):
    '''...

    '''

    def __init__(self, data=None, torsion=None, convention='Berger-Fresse'):
        '''...'''
        # check input data: dict with tuple of int keys
        if data:
            if not (isinstance(data, dict)
                    and all(isinstance(surj, tuple) for surj in data.keys())
                    and all(isinstance(i, int) for i in
                            chain.from_iterable(data.keys()))
                    ):
                raise TypeError(
                    'data type must be dict with tuple of int keys')

        if convention not in {'Berger-Fresse', 'McClure-Smith'}:
            raise TypeError(
                'convention must be either Berger-Fresse or McClure-Smith')

        # initialize element
        self.convention = convention

        super(Surjection_element, self).__init__(data=data, torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def zero(self):
        '''...'''
        return type(self)(torsion=self.torsion,
                          convention=self.convention)

    def create(self, other):
        '''...'''
        return type(self)(other, torsion=self.torsion,
                          convention=self.convention)

    @property
    def arity(self):
        '''...'''
        if not self:
            return None
        arities = set(max(surj) for surj in self.keys())
        if len(arities) == 1:
            return arities.pop()
        return None

    @property
    def degree(self):
        '''...'''
        if not self:
            return None
        degs = set(len(surj) - max(surj) for surj in self.keys())
        if len(degs) == 1:
            return degs.pop()
        return None

    @property
    def complexity(self):
        '''returns the complexity of an element in the Surjection operad

        >>> Surjection_element({(1, 2, 1, 3, 1, 3, 2, 3): 1}).complexity
        3
        >>> Surjection_element({(1, 2, 1):1, (1, 2): 1}).complexity
        2

        '''
        complexities = [0]
        for surjection in self.keys():
            for i, j in combinations(range(1, max(surjection) + 1), 2):
                r = tuple(k for k in surjection if k == i or k == j)
                cpxty = len([p for p, q in pairwise(r) if p != q])
                complexities.append(cpxty)

        return max(complexities)

    def boundary(self):
        '''boundary of self

        >>> s = Surjection_element({(1,2,1,3,1,3,2,3): 1}, torsion=2)
        >>> print(s.boundary())
        (2,1,3,1,3,2,3) + (1,2,3,1,3,2,3) + (1,2,1,3,1,2,3) + (1,2,1,3,1,3,2)
        >>> print(s.boundary().boundary())
        0

        >>> s = Surjection_element({(1,2,1,3,1,3,2,3): 1})
        >>> print(s.boundary())
        (2,1,3,1,3,2,3) + (1,2,3,1,3,2,3) + (1,2,1,3,1,2,3) - (1,2,1,3,1,3,2)
        >>> print(s.boundary().boundary())
        0

        >>> s = Surjection_element({(1,2,1,3,1,3,2,3): 1},
        ...                        convention='McClure-Smith')
        >>> print(s.boundary())
        (2,1,3,1,3,2,3) - (1,2,3,1,3,2,3) + (1,2,1,3,1,2,3) - (1,2,1,3,1,3,2)
        >>> print(s.boundary().boundary())
        0

        '''
        answer = self.zero()

        if self.torsion == 2:
            for k in self.keys():
                for idx in range(0, len(k)):
                    bdry_summand = k[:idx] + k[idx + 1:]
                    if k[idx] in bdry_summand:
                        answer += self.create({bdry_summand: 1})
            return answer

        if self.convention == 'Berger-Fresse':
            for k, v in self.items():
                # determining the signs of the summands
                signs = {}
                alternating_sign = 1
                for idx, i in enumerate(k):
                    if i in k[idx + 1:]:
                        signs[idx] = alternating_sign
                        alternating_sign *= (-1)
                    elif i in k[:idx]:
                        occurs = (pos for pos, j in enumerate(k[:idx]) if i == j)
                        signs[idx] = signs[max(occurs)] * (-1)
                    else:
                        signs[idx] = 0

                # computing the summands
                for idx in range(0, len(k)):
                    bdry_summand = k[:idx] + k[idx + 1:]
                    if k[idx] in bdry_summand:
                        answer += self.create({bdry_summand: signs[idx] * v})

        if self.convention == 'McClure-Smith':
            for k, v in self.items():
                sign = 1
                for i in range(1, max(k) + 1):
                    for idx in (idx for idx, j in enumerate(k) if j == i):
                        new_k = k[:idx] + k[idx + 1:]
                        if k[idx] in new_k:
                            answer += answer.create({new_k: v * sign})
                        sign *= -1
                    sign *= -1

        return answer

    def __rmul__(self, other):
        '''Left action by the appropriate symmetric group ring.

        # chain map checks:

        >>> rho = SymmetricModule.rotation_element(3)
        >>> surj = Surjection_element({(1, 2, 3, 1, 2): 1}, \
                                       convention='Berger-Fresse')
        >>> x, y = (rho * surj).boundary(), rho * surj.boundary()
        >>> x == y
        True

        >>> rho = SymmetricModule.rotation_element(3)
        >>> surj = Surjection_element({(1, 2, 3, 1, 3): 1}, \
                                       convention='McClure-Smith')
        >>> x, y = (rho * surj).boundary(), rho * surj.boundary()
        >>> x == y
        True

        '''
        def sign(perm, surj, convention):
            if convention == 'Berger-Fresse':
                return 1
            assert convention == 'McClure-Smith'
            signs = {0: 1, 1: -1}
            weights = [surj.count(i) - 1 for
                       i in range(1, max(surj) + 1)]
            answer = 0
            for idx, i in enumerate(perm):
                right = [weights[perm.index(j)] for
                         j in perm[idx + 1:] if i > j]
                answer += sum(right) * weights[idx]
            return signs[answer % 2]

        if isinstance(other, int):
            return super().__rmul__(other)

        if isinstance(other, SymmetricGroup_element):
            return SymmetricModule_element({other: 1}, torsion=self.torsion)

        if not isinstance(other, SymmetricModule_element):
            raise TypeError(f'right mult. by type int or \
                SymmetricModule_element not {type(other)}')

        if self.torsion != other.torsion:
            raise TorsionError

        if self.arity != other.arity:
            raise ArityError

        answer = self.zero()
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            new_key = tuple(k2[i - 1] for i in k1)
            new_sign = sign(k2, k1, self.convention)
            answer += self.create({new_key: new_sign * v1 * v2})
        return answer

    def orbit(self, representation='trivial'):
        ''' Returns the preferred element in the symmetric orbit of an element

        >>> s = Surjection_element({(1, 3, 2): 1})
        >>> print(s.orbit(representation='trivial'))
        (1,2,3)
        >>> print(s.orbit(representation='sign'))
        - (1,2,3)

        >>> s = Surjection_element({(2, 1, 2, 1): 1}, \
                                    convention='McClure-Smith')
        >>> print(s.orbit())
        - (1,2,1,2)

        '''
        def sign(permutation, representation):
            if representation == 'trivial':
                return 1
            if representation == 'sign':
                return permutation.sign

        answer = self.zero()
        for k, v in self.items():
            seen = []
            for i in k:
                if i not in seen:
                    seen.append(i)
            permutation = SymmetricGroup_element(seen).inverse()
            new_v = sign(permutation, representation) * v
            answer += permutation * self.create({k: new_v})

        return answer

    def __call__(self, other):
        '''Action on an basis element in the normalized chains of a standard
        cube or simplex represented by an arity 1 element in the (cubical)
        Eilenberg-Zilber operad.

        Examples
        --------

        # chain map check

        >>> s = Surjection_element({(3, 2, 1, 3, 1, 2): 1},\
                                       convention='McClure-Smith')
        >>> x = EilenbergZilber.standard_element(3)
        >>> ds_x = s.boundary()(x)
        >>> d_sx = s(x).boundary()
        >>> sdx = s(x.boundary())
        >>> d_sx - ((-1)**(s.degree)) * sdx == ds_x
        True

        >>> y = EilenbergZilber.standard_element(3)
        >>> ds_y = s.boundary()(y)
        >>> d_sy = s(y).boundary()
        >>> sdy = s(y.boundary())
        >>> d_sy - ((-1)**(s.degree)) * sdy == ds_y
        True

        '''
        def _sign(k1, k2):
            '''...
            '''
            def ordering_sign(permu, weights):
                '''Returns the exponent of the Koszul sign of the given
                permutation acting on the elements of degrees given by the
                list of weights

                '''
                sign_exp = 0
                for idx, j in enumerate(permu):
                    to_add = [weights[permu.index(i)] for
                              i in permu[idx + 1:] if i < j]
                    sign_exp += weights[idx] * sum(to_add)
                return sign_exp % 2

            def action_sign(ordered_k1, ordered_weights):
                '''Given a ordered tuple [1,..,1, 2,...,2, ..., r,...,r]
                and weights [w_1, w_2, ..., w_{r+d}] of the same length, gives
                the kozul sign obtained by inserting from the left a weight 1
                operator between equal consecutive elements.

                '''
                sign_exp = 0
                for idx, (i, j) in enumerate(pairwise(ordered_k1)):
                    if i == j:
                        sign_exp += sum(ordered_weights[:idx + 1])
                return sign_exp % 2

            sign_exp = 0
            weights = [e.dimension % 2 for e in k2]
            inv_ordering_permu = [pair[0] for pair in
                                  sorted(enumerate(k1), key=itemgetter(1))]
            ordering_permu = tuple(inv_ordering_permu.index(i)
                                   for i in range(len(inv_ordering_permu)))
            sign_exp += ordering_sign(ordering_permu, weights)
            ordered_k1 = list(sorted(k1))
            ordered_weights = [weights[i] for i in inv_ordering_permu]
            sign_exp += action_sign(ordered_k1, ordered_weights)
            return (-1)**sign_exp

        def _simplicial(self, other):
            '''...'''
            answer = other.zero()
            pre_join = other.iterated_diagonal(self.arity + self.degree - 1)
            for (k1, v1), (k2, v2) in product(self.items(), pre_join.items()):
                new_k = []
                zero_summand = False
                for i in range(1, max(k1) + 1):
                    to_join = (spx for idx, spx in enumerate(k2)
                               if k1[idx] == i)
                    joined = Simplex(reduce(lambda x, y: x + y, to_join))
                    if joined.is_degenerate():
                        zero_summand = True
                        break
                    new_k.append(joined)

                if not zero_summand:
                    if self.torsion == 2:
                        sign = 1
                    else:
                        sign = _sign(k1, k2)
                    answer += answer.create({tuple(new_k): sign * v1 * v2})
            return answer

        def _cubical(self, other):
            '''...'''
            answer = other.zero()
            pre_join = other.iterated_diagonal(self.arity + self.degree - 1)
            for (k1, v1), (k2, v2) in product(self.items(), pre_join.items()):
                to_dist = []
                zero_summand = False
                for i in range(1, max(k1) + 1):
                    key_to_join = tuple(cube for idx, cube in enumerate(k2)
                                        if k1[idx] == i)
                    joined = other.create({key_to_join: 1}).join()
                    if not joined:
                        zero_summand = True
                        break
                    to_dist.append(joined)

                if not zero_summand:
                    if self.torsion == 2:
                        sign = 1
                    else:
                        sign = _sign(k1, k2)

                    items_to_dist = [summand.items() for summand in to_dist]
                    for pairs in product(*items_to_dist):
                        new_k = reduce(lambda x, y: x + y, (pair[0] for pair in pairs))
                        new_v = reduce(lambda x, y: x * y, (pair[1] for pair in pairs))
                        to_add = answer.create({tuple(new_k): sign * new_v * v1 * v2})
                        answer += to_add

            return answer

        if not self or not other:
            return other.zero()

        if other.arity != 1:
            raise TypeError(f'action only on arity 1, not {other.arity}')

        if self.degree is None or self.arity is None:
            raise TypeError('defined for homogeneous surjections')

        if self.torsion != other.torsion:
            raise TorsionError

        if isinstance(other, EilenbergZilber_element):
            if self.convention != 'McClure-Smith':
                raise NotImplementedError
            return _simplicial(self, other)

        elif isinstance(other, CubicalEilenbergZilber_element):
            return _cubical(self, other)

        else:
            raise NotImplementedError

    def _reduce_rep(self):
        '''Sets to 0 all degenerate surjections.'''
        # remove non-surjections
        zeros = list()
        for k in self.keys():
            if set(k) != set(range(1, max(k) + 1)):
                zeros.append(k)
        for k in zeros:
            del self[k]

        # removes keys w/ equal consecutive values
        for k, v in self.items():
            for i in range(len(k) - 1):
                if k[i] == k[i + 1]:
                    self[k] = 0

        super()._reduce_rep()


class Surjection():
    '''Class producing Surjection elements of special interest.'''

    @staticmethod
    def steenrod_product(arity, degree, torsion=None,
                         convention='Berger-Fresse'):
        '''Returns a surjection element representing the Steenrod
        product in the given arity and degree.

        Constructed recursively by mapping the minimal resolution W(r)
        of Z[S_r] to Surj(r). We use the chain homotopy equivalence
        of Surj(r) and Z defined using the chain contraction (i, p, s)
        relating Surj(r-1) and Surj(r).

        Parameters
        ----------
        arity : int
        Arity of the complex considered, Surj(arity).

        degree : int
        degree of the element considered Surj(arity)_degree.

        Examples
        --------

        # chain map checks:

        >>> t = SymmetricModule.transposition_element(6)
        >>> x = Surjection.steenrod_product(6, 3).boundary()
        >>> y = t * Surjection.steenrod_product(6, 2)
        >>> print(x == y)
        True

        >>> n = SymmetricModule.norm_element(5, torsion=7)
        >>> x = Surjection.steenrod_product(5, 6, torsion=7).boundary()
        >>> y = n * Surjection.steenrod_product(5, 5, torsion=7)
        >>> print(x == y)
        True

        >>> t = SymmetricModule.transposition_element(5)
        >>> x = Surjection.steenrod_product(5, 3, \
            convention='McClure-Smith').boundary()
        >>> y = t * Surjection.steenrod_product(5, 2, \
            convention='McClure-Smith')
        >>> print(x == y)
        True


        '''

        def i(surj, iterate=1):
            '''Inclusion of Surj(r) into Surj(r+1) by appending 1 at
            the start of basis elements and raising the value of all
            other entries by 1.'''

            if iterate == 1:
                answer = surj.zero()
                for k, v in surj.items():
                    answer += answer.create(
                        {(1,) + tuple(j + 1 for j in k): v})
                return answer
            if iterate > 1:
                return i(i(surj, iterate=iterate - 1))

        def p(surj, iterate=1):
            '''Projection of Surj(r) to Surj(r-1) by removing 1
            from a basis element with only one occurence of value 1
            and substracting 1 from all other entries.

            '''
            if iterate == 1:
                answer = surj.zero()
                for k, v in surj.items():
                    if k.count(1) == 1:
                        idx = k.index(1)
                        new_k = (tuple(j - 1 for j in k[:idx]) +
                                 tuple(j - 1 for j in k[idx + 1:]))
                        answer += answer.create({new_k: v})
                return answer
            if iterate > 1:
                return p(p(surj, iterate=iterate - 1))

        def s(surj):
            '''Chain homotopy from the identity to the composition pi, i.e.
            id - ip = ds + sd'''

            answer = surj.zero()
            for k, v in surj.items():
                answer += answer.create({(1,) + tuple(j for j in k): v})
            return answer

        def h(surj):
            '''Chain homotopy from the identiy to i...i p..p in Surj(r)
            realizing its contractibility to Surj(1).

            '''
            answer = s(surj)
            for r in range(1, arity - 1):
                answer += i(s(p(surj, r)), r)
            return answer

        operators = {
            0: SymmetricModule.norm_element(arity),
            1: SymmetricModule.transposition_element(arity)
        }

        def psi(arity, degree, convention=convention):
            '''Recursive definition of the steenrod product over the integers.'''

            if degree == 0:
                return Surjection_element({tuple(range(1, arity + 1)): 1},
                                          convention=convention)
            else:
                previous = psi(arity, degree - 1, convention=convention)
                acted_on = operators[degree % 2] * previous
                answer = h(acted_on)
                return answer

        integral_answer = psi(arity, degree, convention=convention)
        if torsion:
            integral_answer.set_torsion(torsion)
        return integral_answer

    def steenrod_operation(p, s, q, bockstein=False):
        '''Models a chain level representative of P_s or bP_s over the prime p
        acting on an element of degree n'''

        # input check
        if not all(isinstance(i, int) for i in {p, s, q}):
            raise TypeError('initialize with three int p,s,n')
        if not isinstance(bockstein, bool):
            raise TypeError('bockstein must be a boolean')
        if p == 2 and bockstein:
            raise TypeError('bP only defined for odd primes')

        if p == 2:
            coeff = 1
            d = s - q
            if d < 0:
                return Surjection_element(torsion=p)

        else:
            b = int(bockstein)
            # Serre convention: v(2j)=(-1)^j & v(2j+1)=v(2j)*m! w/ m=(p-1)/2
            coeff = (-1)**(floor(q / 2) + s)
            if q / 2 - floor(q / 2):
                coeff *= factorial((p - 1) / 2)
            # degree of the element: (2s-q)(p-1)-b
            d = (2 * s - q) * (p - 1) - b
            if d < 0:
                return Surjection_element(torsion=p)

        return int(coeff) * Surjection.steenrod_product(
            p, d, torsion=p, convention='McClure-Smith')

    @staticmethod
    def basis(arity, degree, complexity=None):
        ''' Returns the list of tuples forming the basis of the surjection
        operad in the given degree, arity and complexity

        >>> basis = sorted(Surjection.basis(3, 5, 3))
        >>> for s in basis:
        ...     print(s)
        (1, 2, 1, 3, 1, 3, 2, 3)
        (1, 3, 1, 2, 1, 2, 3, 2)
        (2, 1, 2, 3, 2, 3, 1, 3)
        (2, 3, 2, 1, 2, 1, 3, 1)
        (3, 1, 3, 2, 3, 2, 1, 2)
        (3, 2, 3, 1, 3, 1, 2, 1)
        '''

        if complexity is None:
            complexity = degree + 1
        a, d, c = arity, degree, complexity
        basis = []
        for s in product(range(1, a + 1), repeat=a + d):
            surj = Surjection_element({s: 1})
            if surj and surj.complexity <= c and surj.arity == a:
                basis.append(s)

        return basis
