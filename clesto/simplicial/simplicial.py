from ..module import Module_element
from ..symmetric import SymmetricRing_element
from ..utils import pairwise
from itertools import chain, product, combinations_with_replacement


class Simplex(tuple):
    """A standard simplex."""

    @property
    def dimension(self):
        """The dimension of self

        Defined as the length of the tuple minus one."""

        return len(self) - 1

    def face(self, i):
        """The i-th face of self

        Obtained by removing the i-th entry of the tuple."""

        return Simplex(self[:i] + self[i + 1:])

    def degeneracy(self, i):
        """The i-th degeneracy of self

        Obtained by repeating the i-th entry of the tuple."""

        return Simplex(self[:i + 1] + self[i:])

    def coface(self, i):
        """The i-th coface of self

        Obtained by adding 1 to each j-th entries with j
        greater or equal to i."""

        def d_i(i, j): return j + 1 if j >= i else j

        return tuple(d_i(i, j) for j in self)

    def codegeneracy(self, i):
        """The i-th codegeneracy of self

        Obtained by substracting 1 from each j-th entries with j
        greater than i."""

        def s_i(i, j): return j - 1 if j > i else j

        return tuple(s_i(i, j) for j in self)

    def is_degenerate(self):
        """Returns True if self is degenerate

        A simplex is degenerate if it is empty or if containes equal
        consecutive values."""

        conseq_values = any([i == j for i, j in pairwise(self)])
        empty_simplex = (self.dimension == -1)
        return empty_simplex or conseq_values


class SimplicialEZ_element_element(Module_element):
    """Element in the Eilenberg-Zilber operad

    This operad is the chain complex of natural transformations from
    the functor of normalized chains to its iterated tensor product
    with itself. An element in arity r is identified with a sequence,
    parametrized by a non-negative integer n, of elements in the
    r-tensor product of the chains on the standard n-simplex which is
    in the kernel of all codegeneracy maps.

    See: J. McClure, and J. Smith. "Multivariable cochain operations and
    little n-cubes." Journal of the American Mathematical Society 16.3
    (2003): 681-704.

    This class models one element of the sequence at the given dimension n.

    """

    def __init__(self, data=None, dimension=None, torsion=None):
        """Initialize an instance of SimplicialEZ_element_element

        Create a new, empty SimplicialEZ_element_element object representing 0,
        and, if given, initialize a SimplicialEZ_element_element from a dict with
        tuple of tuple of int keys and int values.

        Note: this corresponds to the element in the sequence representing a
        natural transformation at the given dimension.

        >>> x = SimplicialEZ_element_element({((0,), (0, 1, 2)): 1, \
                                         ((0, 1), (1, 2)): -1, \
                                         ((0, 1, 2), (2,)): 1})
        """

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

        super(SimplicialEZ_element_element, self).__init__(data=data,
                                                           torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def _latex_(self):
        """Representation in LaTex.

        >>> x = SimplicialEZ_element_element({((0,), (0, 1, 2)): 1, \
                                         ((0, 1), (1, 2)): -1, \
                                         ((0, 1, 2), (2,)): 1})
        >>> print(x._latex_())
        [0] \otimes [0,1,2] - [0,1] \otimes [1,2] + [0,1,2] \otimes [2]

        """
        string = str(self)
        string = string.replace(',),(', '] \otimes [')
        string = string.replace('),(', '] \otimes [')
        string = string.replace('((', '[')
        string = string.replace(',))', ']').replace('),)', ']')
        string = string.replace('))', ']')

        return string

    @property
    def arity(self):
        """Arity of self

        Defined as None if self is not homogeneous. The arity of a basis
        element corresponds to the number of simplices it contains.

        >>> x = SimplicialEZ_element_element({((0,), (0, 1, 2)): 1})
        >>> x.arity
        2

        """
        arities = set(len(multispx) for multispx in self.keys())
        if len(arities) != 1:
            return None
        return arities.pop()

    @property
    def degree(self):
        """Degree of self

        Defined as None if self is not homogeneous. The degree of a basis
        element agrees with the sum of the dimension of the simplices it
        contains.

        >>> x = SimplicialEZ_element_element({((0,), (0, 1, 2)): 1})
        >>> x.degree
        2

        """
        degs = {sum(spx.dimension for spx in k) for k in self.keys()}
        if len(degs) != 1:
            return None
        return degs.pop()

    def boundary(self):
        """Boundary of self

        Defined as the boundary of a tensor product of chains complexes.

        >>> x = SimplicialEZ_element_element({((0, 1), (1, 2)): 1})
        >>> print(x.boundary())
        ((1,),(1,2)) - ((0,),(1,2)) - ((0,1),(2,)) + ((0,1),(1,))

        """
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
        """Left action: other * self

        Left multiplication by a symmetric ring element or an integer.

        >>> x = SimplicialEZ_element_element({((0, 1), (1, 2)): 1})
        >>> t = SymmetricRing_element({(2, 1): 1})
        >>> print(t * x)
        - ((1,2),(0,1))
        >>> print(3 * x)
        3((0,1),(1,2))

        """

        def check_input(self, other):
            """Symmetric ring element with same attributes."""
            if not isinstance(other, SymmetricRing_element):
                raise TypeError(f'right mult. by type int or \
                    SymmetricRing_element not {type(other)}')
            if self.torsion != other.torsion:
                raise TypeError('Unequal torsion attribute')
            if self.arity != other.arity:
                raise TypeError('Unequal arity attribute')

        def sign(perm, multispx):
            weights = [spx.dimension % 2 for spx in k1]
            sign_exp = 0
            for idx, i in enumerate(perm):
                right = [weights[perm.index(j)] for
                         j in perm[idx + 1:] if i > j]
                sign_exp += sum(right) * weights[idx]
            return (-1)**(sign_exp % 2)

        if isinstance(other, int):
            return super().__rmul__(other)
        check_input(self, other)
        answer = self.zero()
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            new_key = [None] * len(k2)
            for idx, i in enumerate(k2):
                new_key[i - 1] = k1[idx]
            new_sign = sign(k2, k1)
            answer += self.create({tuple(new_key): new_sign * v1 * v2})
        return answer

    def coface(self, i):
        """Covariant action of the i-th coface."""

        if i > self.dimension:
            raise TypeError('coface out of range')

        answer = self.zero()
        for k, v in self.items():
            new_k = tuple(tuple(spx.coface(i) for spx in k))
            answer += answer.create({new_k: v})
        return answer

    def codegeneracy(self, i):
        """Covariant action of the i-th codegeneracy. It is the zero map
        on any element in the Eilenberg-Zilber operad.

        """
        if i > self.dimension - 1:
            raise TypeError('codegeneracy out of range')

        answer = self.zero()
        for k, v in self.items():
            new_k = tuple(spx.codegeneracy(i) for spx in k)
            answer += answer.create({new_k: v})
        return answer

    def _reduce_rep(self):
        """Sets to 0 the summands with degenerate simplices."""

        for k, v in self.items():
            if any([spx.is_degenerate() for spx in k]):
                self[k] = 0

        super()._reduce_rep()

    def iterated_diagonal(self, times=1, coord=1):
        """Alexander-Whitney chain approximation to the diagonal applied
        n-times.

        Examples
        --------

        # chain map check:

        >>> x = SimplicialEZ_element_element({((0, 1, 2), ): 1})
        >>> dx = x.boundary()
        >>> dx.iterated_diagonal(3) == x.iterated_diagonal(3).boundary()
        True

        """

        if self.degree is None:
            raise TypeError(f'only for homogeneous elements')

        if self.arity < coord:
            raise TypeError(f'arity = {self.arity} < coord = {coord}')

        answer = self.zero()
        for k, v in self.items():
            left, spx, right = k[:coord - 1], k[coord - 1], k[coord:]
            for p in combinations_with_replacement(range(self.degree + 1), times):
                p = (0,) + p + (self.degree,)
                new_k = []
                for i, j in pairwise(p):
                    new_k.append(Simplex(spx[i:j + 1]))
                answer += self.create({left + tuple(new_k) + right: v})
        return answer


class SimplicialEZ_element():
    """Class producing Eilenberg-Zilber elements of special interest."""

    def standard_element(n, torsion=None):
        return SimplicialEZ_element_element({(tuple(range(n + 1)), ): 1},
                                            torsion=torsion)

    def boundary_element(n, torsion=None):
        """..."""

        sign = {0: 1, 1: -1}
        answer = SimplicialEZ_element_element(dimension=n, torsion=torsion)
        spx = Simplex(range(n + 1))
        for i in range(n + 1):
            new_k = (spx[:i] + spx[i + 1:], )
            answer += answer.create({new_k: sign[i % 2]})
        return answer

    def coproduct_element(n, torsion=None):
        """..."""

        answer = SimplicialEZ_element_element(dimension=n, torsion=torsion)
        spx = Simplex(range(n + 1))
        for i in range(1, n + 2):
            new_k = (spx[:i], spx[i - 1:])
            answer += answer.create({new_k: 1})
        return answer
