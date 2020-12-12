from ..basics import Module_element
from ..basics import SymmetricGroup_element
from ..basics import SymmetricRing_element, SymmetricRing

from ..surjection import Surjection_element
from ..utils import partitions
from itertools import chain, product


class BarrattEccles_element(Module_element):
    """Elements in the Barratt-Eccles operad.

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
            return

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

        >>> x = BarrattEccles_element({((1,3,2), (2,3,1)): 1})
        >>> print(x.boundary().boundary())
        0

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

        The preferred representative in the orbit of basis element is one whose first
        symmetric group element if the identity.

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

    def compose(self, *others):
        """WRONG, NOT DONE BY ME. NOT A CHAIN MAP"""

        def _paths(p, q):
            """returns as a list all increasing paths from (0,0) to (p,q)."""

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
            """..."""
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
        """Table reduction of self

        As defined by Berger-Fresse.

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

    def alexander_whitney(self, r=1):
        """..."""

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
        """Returns representative with only non-zero terms.

        """
        for simplex, v in self.items():
            for i in range(len(simplex) - 1):
                if simplex[i] == simplex[i + 1]:
                    self[simplex] = 0

        super()._reduce_rep()


class BarrattEccles():
    """Class producing Barratt-Eccles elements of special interest."""

    @staticmethod
    def steenrod_product(arity, degree, torsion=None):
        """Returns a representative of the requesed Steenrod product

        Constructed recursively by mapping the minimal resolution W(r)
        of Z[S_r] to E(r). We use the chain homotopy equivalence
        of Surj(r) and Z defined using the chain contraction (i, p, s)
        relating Surj(r-1) and Surj(r).

        """
        operators = {
            0: SymmetricRing.norm_element(arity),
            1: SymmetricRing.transposition_element(arity)
        }

        def _psi(arity, degree):
            """Recursive definition of W --> E over the integers."""

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
