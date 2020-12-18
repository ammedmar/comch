from ..module import Module_element

from ..symmetric import SymmetricGroup_element
from ..symmetric import SymmetricRing_element, SymmetricRing

from ..surjection import Surjection_element
from ..utils import partitions, pairwise
from itertools import chain, product


class BarrattEccles_element(Module_element):
    """Elements in the Barratt-Eccles operad

    As defined in:

    [BF]: C. Berger, and B. Fresse. "Combinatorial operad actions on cochains."
    Mathematical Proceedings of the Cambridge Philosophical Society. Vol. 137.
    No. 1. Cambridge University Press, 2004.

    """

    def __init__(self, data=None, torsion=None):
        """Initialize an instance of BarrattEccles_element

        Create a new, empty BarrattEccles_element object representing 0, and,
        if given, initialize a BarrattEccles_element from a dict with tuple of
        tuple of int keys and int values.

        """
        def check_input_data(data):
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

        def prepare_data(data):
            """transform tuples to symmetric group elements."""
            new_data = {}
            for k, v in data.items():
                new_key = tuple(SymmetricGroup_element(pi) for pi in k)
                new_data[new_key] = v
            return new_data

        if data:
            check_input_data(data)
            data = prepare_data(data)

        # initializing element
        super(BarrattEccles_element, self).__init__(data=data,
                                                    torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    @property
    def arity(self):
        """Arity of self

        Defined as None if self is not homogeneous. The arity of a basis
        element agrees with arity of any of the symmetric group elements

        >>> x = BarrattEccles_element({((1,3,2), (2,3,1)): 1})
        >>> x.arity
        3

        """
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
        """Degree of self

        Defined as None if self is not homogeneous. The degree of a basis
        surjection agrees with the cardinality of the tuple minus one.

        >>> x = BarrattEccles_element({((1,3,2), (2,3,1)): 1})
        >>> x.degree
        1

        """
        if not self:
            return None
        degrees = []
        for k in self.keys():
            degrees.append(len(k) - 1)
        if len(set(degrees)) > 1:
            return None
        return degrees.pop()

    @property
    def complexity(self):
        raise NotImplementedError

    def boundary(self):
        """Boundary of self.

        >>> x = BarrattEccles_element({((1,3,2), (2,3,1), (1,2,3)): 1})
        >>> print(x.boundary())
        ((2,3,1),(1,2,3)) - ((1,3,2),(1,2,3)) + ((1,3,2),(2,3,1))

        """
        sign = {0: 1, 1: -1}
        answer = self.zero()
        for spx, coeff in self.items():
            for i in range(len(spx)):
                i_term = {tuple(spx[: i] + spx[i + 1:]): sign[i % 2] * coeff}
                to_add = answer.create(i_term)
                answer += to_add
        answer._reduce_rep()

        return answer

    def __rmul__(self, other):
        """Left action: other * self

        Left multiplication by a symmetric group element or an integer.

        >>> x = BarrattEccles_element({((1,3,2), (2,3,1)): 1})
        >>> print(-x)
        - ((1,3,2),(2,3,1))
        >>> rho = SymmetricRing_element({(2,3,1): 1})
        >>> print(rho * x)
        ((2,1,3),(3,1,2))

        """
        if isinstance(other, int):
            return super().__rmul__(other)

        if not isinstance(other, SymmetricRing_element):
            raise NotImplementedError

        if self.torsion != other.torsion:
            raise TypeError('Unequal torsion attribute')

        if self.arity != other.arity:
            raise TypeError('Unequal arity attribute')

        answer = self.zero()
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            new_k = tuple(k2 * pi for pi in k1)
            answer += self.create({new_k: v1 * v2})
        return answer

    def orbit(self, representation='trivial'):
        """Returns the preferred representative in the orbit of self

        The preferred representative in the orbit of a basis element is one
        whose first symmetric group element if the identity.

        The representation used can be either 'trivial' or 'sign'.

        >>> x = BarrattEccles_element({((1,3,2), (1,2,3)): 1})
        >>> print(x.orbit())
        ((1,2,3),(1,3,2))
        >>> print(x.orbit('sign'))
        - ((1,2,3),(1,3,2))

        """
        if not self:
            return self

        answer = BarrattEccles_element(torsion=self.torsion)
        for k, v in self.items():
            inverse = tuple(k[0].index(i + 1) + 1 for i in range(len(k[0])))
            permutation = SymmetricRing_element({inverse: 1}, torsion=self.torsion)
            if representation == 'sign':
                permutation = k[0].sign * permutation
            answer += permutation * BarrattEccles_element({k: v}, torsion=self.torsion)

        return answer

    def compose(self, other, position):
        """Operadic compositions: self o_position other

        We think of other being inserted into self the pair is ordered:
        self tensor other.

        >>> x = BarrattEccles_element({((1, 2), (2, 1)): 1})
        >>> print(x.compose(x, 1))
        - ((1,2,3),(2,1,3),(3,2,1)) + ((1,2,3),(3,1,2),(3,2,1))

        """

        def check_input(self, other, position):
            """Homogeneous, equal torsion, and position less than arity."""
            if self.torsion != other.torsion:
                raise TypeError('Unequal torsion attribute')
            if not self or not other:
                return self.zero()
            if None in [self.arity, other.arity, self.degree, other.degree]:
                raise TypeError('Defined for homogeneous elements only')
            if self.arity < position:
                raise TypeError(f'Arity {self.arity} < position = {position}')

        def paths(p, q):
            """All increasing paths from (0,0) to (p,q)."""
            if (p, q) == (0, 0):
                return [((0, 0),)]
            answer = list()
            if p > 0:
                west = paths(p - 1, q)
                for path in west:
                    answer.append(path + ((p, q),))
            if q > 0:
                south = paths(p, q - 1)
                for path in south:
                    answer.append(path + ((p, q),))
            return answer

        def sign_of_path(path):
            """Sign of shuffle placing horizontal before vertical lines."""
            vectors = [(a[0] - b[0], a[1] - b[1]) for b, a in pairwise(path)]
            sign_exp = 0
            for idx, vector in enumerate(vectors):
                if vector == (0, 1):
                    sign_exp += len([v for v in vectors[idx + 1:] if v == (1, 0)])
            return (-1) ** (sign_exp)

        check_input(self, other, position)
        answer = self.zero()
        p, q = self.degree, other.degree
        all_paths = paths(p, q)
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            for path in all_paths:
                new_k = []
                for i, j in path:
                    new_k.append(k1[i].compose(k2[j], position))
                sign = 1 if self.torsion == 2 else sign_of_path(path)
                answer += self.create({tuple(new_k): sign * v1 * v2})
        return answer

    def table_reduction(self):
        """Table reduction of self

        >>> b = BarrattEccles_element({((1,2,3,4), (1,4,3,2), (1,2,4,3)): 1})
        >>> print(b.table_reduction())
        (1,2,4,2,4,3) + (1,2,4,3,2,3)

        """
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

    def diagonal(self, r=1):
        """Alexander Whitney diagonal

        Defined on basis elements by
        sum_i [pi_0,...pi_i] tensor [pi_i,...,pi_d]

        """

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

    def _reduce_rep(self):
        """Removes degenerate elements

        Those with consecutive equal permutations or equal to an empty
        permutation.

        """
        for k, v in self.items():
            if not k:
                self[k] = 0
            for per1, per2 in pairwise(k):
                if per1 == per2:
                    self[k] = 0
        super()._reduce_rep()


class BarrattEccles():
    """Class producing Barratt-Eccles elements of special interest."""

    @staticmethod
    def steenrod_product(arity, degree, torsion=None):
        """Returns a representative of the requested Steenrod product

        Constructed recursively by mapping the minimal resolution W(r)
        of Z[S_r] to E(r). We use the chain homotopy equivalence
        of Surj(r) and Z defined using the chain contraction (i, p, s)
        relating Surj(r-1) and Surj(r).

        """
        operators = {
            0: SymmetricRing.norm_element(arity),
            1: SymmetricRing.transposition_element(arity)
        }

        def psi(arity, degree):
            """Recursive definition of W --> E over the integers."""

            if degree == 0:
                return BarrattEccles_element(
                    {(tuple(range(1, arity + 1)),): 1})
            else:
                previous = psi(arity, degree - 1)
                acted_on = operators[degree % 2] * previous
                identity = tuple(range(1, arity + 1))
                new_data = {(identity,) + k: v for k, v in acted_on.items()}
                return BarrattEccles_element(new_data)

        integral_answer = psi(arity, degree)
        if torsion:
            integral_answer.set_torsion(torsion)
        return integral_answer
