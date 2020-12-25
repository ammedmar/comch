from ..module import ModuleElement

from ..symmetric import SymmetricGroupElement
from ..symmetric import SymmetricRingElement, SymmetricRing

from ..surjection import SurjectionElement
from ..utils import partitions, pairwise
from itertools import product


class BarrattEcclesElement(ModuleElement):
    r"""Element in the Barratt-Eccles operad

    For a non-negative integer :math:`r` define the simplicial set
    :math:`E(\mathrm S_r)` by

    .. math::

       \begin{aligned}
       E(\mathrm S_r)_n &=
       \{ (\sigma_0, \dots, \sigma_n)\ |\ \sigma_i \in \mathrm{S}_r\}, \\
       d_i(\sigma_0, \dots, \sigma_n) &=
       (\sigma_0, \dots, \widehat{\sigma}_i, \dots, \sigma_n), \\
       s_i(\sigma_0, \dots, \sigma_n) &=
       (\sigma_0, \dots, \sigma_i, \sigma_i, \dots, \sigma_n),
       \end{aligned}

    corresponding to the unreduced bar construction on the monoid
    :math:`\mathrm S_r`. It is equipped with a left
    :math:`\mathrm S_r`-action defined on basis elements by

    .. math::
        \sigma (\sigma_0, \dots, \sigma_n) =
        (\sigma \sigma_0, \dots, \sigma \sigma_n).

    The chain complex resulting from applying the functor of integral
    normalized chains to it is denoted :math:`\mathcal E(r)`, which
    corresponds to the arity :math:`r` part of the Barratt-Eccles operad.

    PARAMETERS
    ----------

    data : ``dict`` or ``None``, default: ``None``
        Dictionary representing a linear cobination of basis elements.
        Items in the dict correspond with pairs `basis_element: coefficient`.
        Each basis_element must create a ``tuple`` of `SymmetricGroupElement`
        and `coefficient` must be an ``int``.
    torsion : ``int`` or ``string`` 'free', default 'free'
        The torsion of the underlying ring.

    ATTRIBUTES
    ----------

    torsion : int or 'free'
        The torsion of the underlying ring.

    EXAMPLES
    --------

    >>> x = BarrattEcclesElement()
    >>> print(x)
    0
    >>> y = BarrattEcclesElement({((1,3,2), (2,1,3)): -1})
    >>> print(y)
    - ((1,3,2),(2,1,3))

    """

    def __init__(self, data=None, torsion=None):
        if data:
            new_data = {}
            for k, v in data.items():
                new_key = tuple(SymmetricGroupElement(pi) for pi in k)
                new_data[new_key] = v
            data = new_data

        super(BarrattEcclesElement, self).__init__(data=data, torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    @property
    def arity(self):
        """Arity of self

        Defined as None if self is not homogeneous. The arity of a basis
        element agrees with arity of any of the symmetric group elements

        EXAMPLE
        -------

        >>> x = BarrattEcclesElement({((1,3,2), (2,3,1)): 1})
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

        EXAMPLE
        -------

        >>> x = BarrattEcclesElement({((1,3,2), (2,3,1)): 1})
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

        >>> x = BarrattEcclesElement({((1,3,2), (2,3,1), (1,2,3)): 1})
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
        answer.preferred_rep()

        return answer

    def __rmul__(self, other):
        """Left action: other * self

        Left multiplication by a symmetric group element or an integer.

        >>> x = BarrattEcclesElement({((1,3,2), (2,3,1)): 1})
        >>> print(-x)
        - ((1,3,2),(2,3,1))
        >>> rho = SymmetricRingElement({(2,3,1): 1})
        >>> print(rho * x)
        ((2,1,3),(3,1,2))

        """
        if isinstance(other, int):
            return super().__rmul__(other)

        if not isinstance(other, SymmetricRingElement):
            raise NotImplementedError

        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')

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

        >>> x = BarrattEcclesElement({((1,3,2), (1,2,3)): 1})
        >>> print(x.orbit())
        ((1,2,3),(1,3,2))
        >>> print(x.orbit('sign'))
        - ((1,2,3),(1,3,2))

        """
        if not self:
            return self

        answer = BarrattEcclesElement(torsion=self.torsion)
        for k, v in self.items():
            inverse = tuple(k[0].index(i + 1) + 1 for i in range(len(k[0])))
            permutation = SymmetricRingElement({inverse: 1}, torsion=self.torsion)
            if representation == 'sign':
                permutation = k[0].sign * permutation
            answer += permutation * BarrattEcclesElement({k: v}, torsion=self.torsion)

        return answer

    def compose(self, other, position):
        """Operadic compositions: self o_position other

        We think of other being inserted into self the pair is ordered:
        self tensor other.

        >>> x = BarrattEcclesElement({((1, 2), (2, 1)): 1})
        >>> print(x.compose(x, 1))
        - ((1,2,3),(2,1,3),(3,2,1)) + ((1,2,3),(3,1,2),(3,2,1))

        """

        def check_input(self, other, position):
            """Homogeneous, equal torsion, and position less than arity."""
            if self.torsion != other.torsion:
                raise TypeError('only defined for equal attribute torsion')
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

        >>> b = BarrattEcclesElement({((1,2,3,4), (1,4,3,2), (1,2,4,3)): 1})
        >>> print(b.table_reduction())
        (1,2,4,2,4,3) + (1,2,4,3,2,3)

        """
        answer = SurjectionElement(torsion=self.torsion,
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

        answer.preferred_rep()

        return answer

    def diagonal(self, r=1):
        """Alexander Whitney diagonal

        Defined on basis elements by
        sum_i [pi_0,...pi_i] tensor [pi_i,...,pi_d]

        """

        def split(multispx):
            a, b = multispx[0], multispx[1:]
            return set((a[:i + 1], a[i:]) + b for i in range(len(a)))

        answer = ModuleElement(torsion=self.torsion)
        for k, v in self.items():
            to_add = {(k,)}
            for s in range(1, r + 1):
                to_add = set.union(*(split(multispx) for multispx in to_add))
            answer += ModuleElement(
                {multispx: v for multispx in to_add}).copy_attrs_from(self)

        return answer

    def preferred_rep(self):
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
        super().preferred_rep()


class BarrattEccles:
    """Produces Barratt-Eccles elements of interest."""

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
                return BarrattEcclesElement(
                    {(tuple(range(1, arity + 1)),): 1})
            else:
                previous = psi(arity, degree - 1)
                acted_on = operators[degree % 2] * previous
                identity = tuple(range(1, arity + 1))
                new_data = {(identity,) + k: v for k, v in acted_on.items()}
                return BarrattEcclesElement(new_data)

        integral_answer = psi(arity, degree)
        if torsion:
            integral_answer.set_torsion(torsion)
        return integral_answer
