from ..free_module import FreeModuleElement

from ..symmetric import SymmetricGroupElement
from ..symmetric import SymmetricRingElement, SymmetricRing

from ..surjection import SurjectionElement
from ..simplicial import SimplicialElement
from ..cubical import CubicalElement
from ..utils import partitions, pairwise
from itertools import product, combinations
from math import comb, factorial


class BarrattEcclesElement(FreeModuleElement):
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

    REFERENCES
    ----------
    [BF]: C. Berger, and B. Fresse. "Combinatorial operad actions on cochains."
    Mathematical Proceedings of the Cambridge Philosophical Society. Vol. 137.
    No. 1. Cambridge University Press, 2004.

    """

    def __init__(self, data=None, torsion=None):
        r"""Initializes *self*.

        PARAMETERS
        ----------
        data : :class:`dict` or ``None``, default: ``None``
            Dictionary representing a linear combination of basis elements. Items
            in the dictionary correspond with pairs `basis_element: coefficient`.
            Each basis_element must create a :class:`tuple` of
            :class:`commch.symmetric.SymmetricGroupElement` and `coefficient` must
            be an :class:`int`.
        torsion : :class:`int` or ``None``, default ``None``
            The integer :math:`n` of the ring :math:`\mathbb Z/n \mathbb Z`;
            :math:`n = 0` means :math:`\mathbb Z`, and negative values are
            normalized to their absolute value.

        EXAMPLES
        --------
        >>> x = BarrattEcclesElement()
        >>> print(x)
        0
        >>> y = BarrattEcclesElement({((1,3,2), (2,1,3)): -1})
        >>> print(y)
        - ((1,3,2),(2,1,3))

        """
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
        """Arity of *self*.

        Defined as ``None`` if *self* is not homogeneous. The arity of a basis
        Barratt-Eccles element agrees with arity of any of the symmetric group
        elements contained in it.

        RETURNS
        _______
        :class:`int`
            The arity of *self*.

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
            arities_in_k = set(perm.arity for perm in k)
            arities |= arities_in_k  # in place union
        if len(arities) == 1:
            return arities.pop()
        return None

    @property
    def degree(self):
        """Degree of *self*.

        Defined as ``None`` if *self* is not homogeneous. The degree of a basis
        Barratt-Eccles element agrees with the cardinality of the tuple minus one.

        RETURNS
        _______
        :class:`int`
            The degree of *self*.

        EXAMPLE
        -------
        >>> x = BarrattEcclesElement({((1,3,2), (2,3,1)): 1})
        >>> x.degree
        1

        """
        if not self:
            return None
        degrees = set(len(k) - 1 for k in self.keys())
        if len(set(degrees)) == 1:
            return degrees.pop()
        return None

    @property
    def complexity(self):
        r"""Complexity of *self*.

        Defined as ``None`` if *self* is not homogeneous.
        The complexity of a finite binary sequence of elements in
        :math:`\Sigma_2` is defined as the number of consecutive distinct
        elements in it. For example, :math:`((12),(21),(21),(12))` and
        :math:`((12),(12),(12),(21))` have complexities 2 and 1 respectively.
        For any basis Barratt-Eccles element, and any pair of positive integers
        :math:`i < j` less than its arity, we can form a sequence as above by
        precomposing each permutation by the order-preserving inclusion sending
        :math:`1` and :math:`2` respectively to :math:`i` and :math:`j`. The
        complexity of a basis Barratt-Eccles element is defined as the maximum
        over $i < j$ of the complexities of these. Notice that for arity 2,
        the complexity of a basis element is equal to its degree. It is proven
        in [BF] that the subcomplex generated by basis Barratt-Eccles elements
        of complexity at most :math:`n` define a suboperad of :math:`\mathcal E`
        modeling an :math:`E_{n+1}`-operad.

        RETURNS
        _______
        :class:`int`
            The complexity of *self*.

        EXAMPLE
        -------
        >>> BarrattEcclesElement({((1,2,3), (1,3,2), (1,2,3)): 1}).complexity
        2

        """
        complexity = 0
        for k in self.keys():
            for i, j in combinations(range(1, max(k[0]) + 1), 2):
                seq = (filter(lambda x: x == i or x == j, perm) for perm in k)
                cpxty = len([p for p, q in pairwise(seq) if p != q])
                complexity = max(cpxty, complexity)
        return complexity

    def boundary(self):
        r"""Boundary of *self*.

        It is defined as the alternating sum of the face maps. Explicitly,
        for basis Barrat-Eccles elements, we have

        .. math::

           \begin{equation*}
           \partial(\sigma_0, \dots, \sigma_n) =
           \sum_{i=0}^n(\sigma_0, \dots, \widehat{\sigma}_i, \dots, \sigma_n).
           \end{equation*}

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The boundary of *self*.

        EXAMPLE
        -------
        >>> x = BarrattEcclesElement({((1,3,2), (2,3,1), (1,2,3)): 1})
        >>> print(x.boundary())
        ((2,3,1),(1,2,3)) - ((1,3,2),(1,2,3)) + ((1,3,2),(2,3,1))

        """
        answer = self.zero()
        for spx, coeff in self.items():
            for i in range(len(spx)):
                i_term = {tuple(spx[:i] + spx[i + 1:]): (-1) ** (i % 2) * coeff}
                to_add = answer.create(i_term)
                answer += to_add
        answer.preferred_rep()
        return answer

    def __rmul__(self, other):
        r"""Left action: *other* ``*`` *self*.

        Left multiplication by a symmetric group element or an integer.
        Defined on basis Barratt-Eccles elements by acting coordinatewise.
        Explicitly,

        .. math::
            \sigma (\sigma_0, \dots, \sigma_n) =
            (\sigma \sigma_0, \dots, \sigma \sigma_n).

        PARAMETERS
        ----------
        other : :class:`int` or :class:`comch.symmetric.SymmetricRingElement`.
            The element to left act on *self* with.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The product: *other* ``*`` *self*.

        EXAMPLE
        -------
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
        """The preferred representative of the symmetric orbit of *self*.

        The preferred representative in the orbit of a basis element is one
        whose first symmetric group element is the identity.

        The representation used can be either 'trivial' or 'sign'.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The preferred element in the symmetric orbit of *self*.

        EXAMPLE
        -------
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
        r"""Operadic compositions: *self* :math:`\circ_{position}` *other*.

        The operadic composition can be described in terms of the composition of
        symmetric group elements using the Eilenberg-Zilber map. Let us notice
        that at the level of the simplicial set :math:`E` we have compositions
        induced coordinate-wise

        .. math:: {\circ}_{i}: E(r) \times E(s) \to E(r + s - 1).

        We define the composition of :math:`\mathcal E` by precomposing

        .. math::
           N_\bullet(\circ_i) \colon N_\bullet(E(r) \times E(s))
           \longrightarrow
           N_\bullet(E(r + s - 1)) = \mathcal E(r+s-1)

        with the iterated Eilenberg-Zilber map

        .. math::
           \mathcal E(r) \otimes \mathcal E(s) =
           N_\bullet(E(r)) \otimes N_\bullet(E(s))
           \longrightarrow
           N_\bullet(E(r) \times E(s)).

        PARAMETERS
        ----------
        other : :class:`comch.barratt_eccles.BarrattEcclesElement`
            The element to operad compose *self* with.
        position : :class:`int`
            The value at which the composition occurs.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The operadic composition of *self* and *other*.

        EXAMPLE
        -------
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

    def suspension(self):
        raise NotImplementedError

    def table_reduction(self):
        r"""Table reduction of *self*.

        The table reduction morphisms :math:`TR : \mathcal E \to \mathcal X`
        from the Barratt-Eccles to the surjection operad is a surjective weak
        equivalence of operads introduced in [BF]. For a basis Barratt-Eccles
        element :math:`(\sigma_0, \dots, \sigma_n) \in \mathcal E(r)_n` we have
        that

        .. math:: TR(\sigma_0, \dots, \sigma_n) = \sum_{a} s_{a}

        with surjections

        .. math:: s_{a} : \{1, \dots, n+r \} \to \{1, \dots, r\}

        parametrized by all tuples of positive integers
        :math:`a = (a_0, \dots, a_n)` with :math:`a_0 + \cdots + a_n = n + r`.
        For one such tuple :math:`a` we now describe the surjection :math:`s_a`.
        Define recursively

        .. math:: A_{-1} = 0 \qquad A_i = A_{i-1} + a_{i.}

        For :math:`k \in \{1, \dots, n+r\}` we identify
        :math:`i \in \{1, \dots, n\}` such that :math:`A_{i-1} < k \leq A_{i}`
        and define :math:`s_a(k)` to be the :math:`(k - A_{i-1})`-th element in
        :math:`(\sigma_i(1), \dots, \sigma_i(r))` not in

        .. math::
            \big\{s_a(j)\ |\ j < k \text{ and } j \neq A_0, \dots, A_{i-1}\big\}.

        This operad map preserves the little-cubes filtration: it sends the
        suboperad spanned by basis Barratt-Eccles elements of complexity at
        most :math:`n - 1` to the one spanned by basis surjection elements
        of complexity at most :math:`n`, both modeling an
        :math:`E_n`-operad.

        RETURNS
        _______
        :class:`comch.surjection.SurjectionElement`
            The image of *self* under the table reduction morphism

        EXAMPLE
        -------
        >>> b = BarrattEcclesElement({((1,2,3,4), (1,4,3,2), (1,2,4,3)): 1})
        >>> print(b.table_reduction())
        (1,2,4,2,4,3) + (1,2,4,3,2,3)

        SEE ALSO
        --------
        :attr:`comch.surjection.SurjectionElement.complexity`
        :attr:`comch.barratt_eccles.BarrattEcclesElement.complexity`

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
        r"""Alexander-Whitney diagonal.

        The Alexander-Whitney chain approximation to the diagonal is the chain
        map :math:`\Delta \colon \mathcal E \to \mathcal E \otimes \mathcal E`
        defined on a basis element
        :math:`(\sigma_0, \dots, \sigma_n) \in \mathcal E(r)_n` by the formula

        .. math::
            \Delta (\sigma_0, \dots, \sigma_n) =
            \sum_{i=1}^n (\sigma_0, \dots, \sigma_i) \otimes
            (\sigma_i, \dots, \sigma_n).

        It defines a Hopf structure on the Barratt-Eccles operad.

        RETURNS
        _______
        :class:`comch.free_module.FreeModuleElement`
            The free module element representing the image of *self* under
            the diagonal map. The keys are pairs of elements in
            :class:`comch.barratt_eccles.BarrattEcclesElement` coefficient 1.

        EXAMPLE
        -------
        >>> x = BarrattEcclesElement({((1,2), (2,1)):1})
        >>> print(x.diagonal())
        (((1, 2),), ((1, 2), (2, 1))) + (((1, 2), (2, 1)), ((2, 1),))

        """

        def split(multispx):
            a, b = multispx[0], multispx[1:]
            return set((a[:i + 1], a[i:]) + b for i in range(len(a)))

        answer = FreeModuleElement(torsion=self.torsion)
        for k, v in self.items():
            to_add = {(k,)}
            for s in range(1, r + 1):
                to_add = set.union(*(split(multispx) for multispx in to_add))
            answer += answer.create({multispx: v for multispx in to_add})

        return answer

    def __call__(self, other, coord=1):
        r"""The action ``self(other)``.

        The action of *self* on the tensor factor specified by *coord* of an
        element in a tensor product of normalized chains of a standard
        simplex or of a standard cube.

        On simplicial elements, the action realizes the natural
        :math:`\mathcal E`-coalgebra structure on normalized chains
        constructed in [BF]. It is the composition of the table reduction
        morphism :math:`TR \colon \mathcal E \to \mathcal X` and the interval
        cut action of the surjection operad with the Berger-Fresse sign
        convention.

        On cubical elements, a direct action of the Barratt-Eccles operad
        is used. Its construction is forthcoming.

        PARAMETERS
        ----------
        other : :class:`comch.simplicial.SimplicialElement` or \
:class:`comch.cubical.CubicalElement`
            The element to which apply *self* to.
        coord : :class:`int`, default 1.
            The tensor factor of *other* to apply *self* to.

        RETURNS
        _______
        :class:`comch.simplicial.SimplicialElement` or \
:class:`comch.cubical.CubicalElement`
            The action of *self* on *other* at *coord*.

        REFERENCES
        ----------
        [BF]: C. Berger, and B. Fresse. "Combinatorial operad actions on
        cochains." Mathematical Proceedings of the Cambridge Philosophical
        Society. Vol. 137. No. 1. Cambridge University Press, 2004.

        EXAMPLE
        -------
        >>> from comch.simplicial import Simplicial
        >>> x = BarrattEcclesElement({((1, 2), (2, 1)): 1})
        >>> print(x(Simplicial.standard_element(2)))
        - ((0,1,2),(0,1)) + ((0,2),(0,1,2)) - ((0,1,2),(1,2))

        """
        if isinstance(other, SimplicialElement):
            return self.table_reduction()(other, coord)
        if isinstance(other, CubicalElement):
            return self._cubical_action(other, coord)
        raise NotImplementedError

    def _cubical_action(self, other, coord=1):
        """To be constructed: action on the chains of standard cubes."""
        raise NotImplementedError(
            'the action of the Barratt-Eccles operad on cubical chains '
            'is not yet implemented')

    def preferred_rep(self):
        """Preferred representative of *self*.

        Removes pairs `basis_element: coefficient` which satisfy either of:
        1) the basis element has equal consecutive permutations, 2) it is the
        empty permutation, or 3) the coefficient is 0.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The preferred representative of *self*.

        EXAMPLE
        -------
        >>> print(BarrattEcclesElement({((1,2),(1,2)):1, ():1, ((1,2),):0}))
        0
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
    def may_steenrod_map(arity, degree, torsion=None):
        r"""Image of the map from the minimal resolution of the cyclic group.

        Returns the image under the May-Steenrod structure constructed in
        [KMM] of the preferred basis element, of the given degree, of the
        minimal free resolution of the base ring :math:`R` as an
        :math:`R[\mathrm C_r]`-module.

        PARAMETERS
        ----------
        arity : :class:`int`
            The arity considered.
        degree : :class:`int`
            The degree considered.
        torsion : :class:`int` or ``None``, default ``None``
            The integer :math:`n` of the ring :math:`\mathbb Z/n \mathbb Z`;
            :math:`n = 0` means :math:`\mathbb Z`, and negative values are
            normalized to their absolute value.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The image of the basis element under the May-Steenrod structure.

        REFERENCES
        ----------
        [KMM]: R. M. Kaufmann and A. M. Medina-Mardones. "Cochain level
        May-Steenrod operations." Forum Mathematicum 33 (2021), no. 6,
        1507--1526. https://doi.org/10.1515/forum-2020-0296

        SEE ALSO
        --------
        :class:`comch.steenrod.Steenrod`

        """
        operators = {
            0: SymmetricRing.norm_element(arity),
            1: SymmetricRing.transposition_element(arity)
        }

        def psi(arity, degree):
            if degree == 0:
                return BarrattEcclesElement(
                    {(tuple(range(1, arity + 1)),): 1})
            previous = psi(arity, degree - 1)
            acted_on = operators[degree % 2] * previous
            identity = tuple(range(1, arity + 1))
            new_data = {(identity,) + k: v for k, v in acted_on.items()}
            return BarrattEcclesElement(new_data)

        integral_answer = psi(arity, degree)
        if torsion is not None:
            integral_answer.set_torsion(torsion)
        return integral_answer

    @staticmethod
    def steenrod_operation_representative(p, s, q, bockstein=False):
        r"""Representative of :math:`P_s` or :math:`\beta P_s`.

        Let :math:`A` be an arbitrary algebra over the Barratt-Eccles operad
        :math:`\mathcal X`. For a prime :math:`p`, an
        integer :math:`s`, and a mod-:math:`p` cocycle :math:`a` of degree
        :math:`q`, this method constructs a Barratt-Eccles element whose action on
        :math:`a^{\otimes p}` is a cocycle representing
        :math:`P_s[a]`. When :math:`p` is odd and ``bockstein=True``, it
        instead constructs a representative of :math:`\beta P_s[a]`.

        PARAMETERS
        ----------
        p : :class:`int`
            The prime considered.
        s : :class:`int`
            The subscript of the Steenrod operation.
        q : :class:`int`
            The degree of the class acted on.
        bockstein : :class:`bool`, default ``False``
            Determines the use of the Bockstein homomorphism.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The Barratt-Eccles element representative.

        REFERENCES
        ----------
        [KMM]: R. M. Kaufmann and A. M. Medina-Mardones. "Cochain level
        May-Steenrod operations." Forum Mathematicum 33 (2021), no. 6,
        1507--1526. https://doi.org/10.1515/forum-2020-0296

        SEE ALSO
        --------
        :class:`comch.steenrod.Steenrod`

        """
        if not all(isinstance(i, int) for i in {p, s, q}):
            raise TypeError('initialize with three int p,s,n')
        if not isinstance(bockstein, bool):
            raise TypeError('bockstein must be a boolean')
        if p == 2 and bockstein:
            raise TypeError('bP only defined for odd primes')

        if p == 2:
            coeff = 1
            degree = s - q
            if degree < 0:
                return BarrattEcclesElement(torsion=2)
        else:
            b = int(bockstein)
            coeff = (-1) ** (q // 2 + s)
            if q % 2:
                coeff *= factorial((p - 1) // 2)
            degree = (2 * s - q) * (p - 1) - b
            if degree < 0:
                return BarrattEcclesElement(torsion=p)

        element = BarrattEccles.may_steenrod_map(p, degree, torsion=p)
        return int(coeff) * element

    @staticmethod
    def cartan_homotopy(i):
        r"""Representative of a Cartan :math:`i`-coboundary.

        Let :math:`A` be an algebra over the Barratt-Eccles operad
        :math:`\mathcal E`. For non-negative :math:`i` and mod-:math:`2`
        cocycles :math:`\alpha,\beta \in A`, this method constructs an arity-4
        Barratt-Eccles element whose action on

        .. math::
            \alpha \otimes \alpha \otimes \beta \otimes \beta

        defines a Cartan :math:`i`-coboundary
        :math:`\zeta_i(\alpha \otimes \beta)`. It satisfies

        .. math::
            \begin{aligned}
            \delta\zeta_i(\alpha \otimes \beta)
            ={}&
            (\alpha \smile_0 \beta)
            \smile_i
            (\alpha \smile_0 \beta) \\
            &+
            \sum_{j+k=i}
            (\alpha \smile_j \alpha)
            \smile_0
            (\beta \smile_k \beta).
            \end{aligned}

        Consequently, :math:`\zeta_i` gives an explicit cochain-level witness
        of the Cartan formula for Steenrod squares.

        More precisely, let

        .. math::
            \widetilde{x}_i =
            (e,(12),e,\ldots,(12)^i) \in \mathcal E(2).

        This method returns the normalized Barratt-Eccles element

        .. math::
            (H_1 + H_2)(\widetilde{x}_i) \in \mathcal E(4),

        where :math:`H_1` and :math:`H_2` are the degree-one maps constructed
        in [Med20]. Its action defines the cochain above by

        .. math::
            \zeta_i(\alpha \otimes \beta) =
            (H_1 + H_2)(\widetilde{x}_i)
            (\alpha \otimes \alpha \otimes \beta \otimes \beta).

        PARAMETERS
        ----------
        i : :class:`int`
            A non-negative degree.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The arity-4, degree-:math:`i+1` Cartan homotopy element over
            :math:`\mathbb F_2`.

        REFERENCES
        ----------
        [Med20]: A. M. Medina-Mardones. "An effective proof of the Cartan
        formula: the even prime." Journal of Pure and Applied Algebra 224
        (2020), no. 12, 106444.
        https://doi.org/10.1016/j.jpaa.2020.106444

        EXAMPLE
        -------
        >>> print(BarrattEccles.cartan_homotopy(0))
        ((1,3,2,4),(1,2,3,4))

        """
        if not isinstance(i, int) or isinstance(i, bool):
            raise TypeError('i must be an integer')
        if i < 0:
            raise ValueError('i must be non-negative')

        identity = (1, 2, 3, 4)
        transposition = (2, 1)
        x = tuple((1, 2) if j % 2 == 0 else transposition
                  for j in range(i + 1))

        first_map = {
            (1, 2): (1, 3, 2, 4),
            transposition: (2, 4, 1, 3),
        }
        second_map = {
            (1, 2): identity,
            transposition: (2, 1, 4, 3),
        }
        first_values = tuple(first_map[permutation] for permutation in x)
        second_values = tuple(second_map[permutation] for permutation in x)

        # H1 is the prism homotopy between the conjugate inclusions f and g.
        answer = BarrattEcclesElement(torsion=2)
        for j in range(i + 1):
            simplex = first_values[:j + 1] + second_values[j:]
            answer += answer.create({simplex: 1})

        # H2 applies the closed Shih homotopy formula before composition.
        composition = {
            ((1, 2), (1, 2)): identity,
            ((1, 2), transposition): (1, 2, 4, 3),
            (transposition, (1, 2)): (2, 1, 3, 4),
            (transposition, transposition): (2, 1, 4, 3),
        }

        def face(simplex, position):
            return simplex[:position] + simplex[position + 1:]

        def degeneracy(simplex, position):
            return simplex[:position + 1] + simplex[position:]

        for p in range(i):
            for q in range(i - p):
                m = i - p - q
                for v_indices in combinations(range(p + q + 1), p):
                    w_indices = tuple(j for j in range(p + q + 1)
                                      if j not in v_indices)

                    first = x
                    for j in reversed(range(i - p + 1, i + 1)):
                        first = face(first, j)
                    first = degeneracy(first, m - 1)
                    for j in v_indices:
                        first = degeneracy(first, j + m)

                    second = x
                    for j in reversed(range(i - p - q, i - p)):
                        second = face(second, j)
                    for j in w_indices:
                        second = degeneracy(second, j + m)

                    simplex = tuple(composition[pair]
                                    for pair in zip(first, second))
                    answer += answer.create({simplex: 1})

        return answer

    @staticmethod
    def adem_homotopy(q, p):
        r"""Homotopy component of an Adem coboundary at the even prime.

        This method implements the chain-level homotopy underlying the
        cochain proof of the mod-2 Adem relations in [BMM]. For integers
        :math:`q > p \geq 0`, it returns the Barratt-Eccles element

        .. math::

           J_\Psi(\widetilde{x}_q \times \widetilde{x}_p)
           + J_\Psi(\widetilde{x}_p \times \widetilde{x}_q)
           + J_{(23)}(\widetilde{x}_q \times \widetilde{x}_p)
           \in \mathcal E(4),

        where :math:`\Psi = EZ \circ \Phi \circ AW` is the map of Section 3.4
        of [BMM] and :math:`J_{(23)}` is the homotopy associated to the inner
        automorphism by the transposition :math:`(23)`.

        PARAMETERS
        ----------
        q, p : :class:`int`
            Non-negative integers with ``q > p``.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The arity-4, degree-:math:`p+q+1` Adem homotopy element over
            :math:`\mathbb F_2`.

        REFERENCES
        ----------
        [BMM]: G. Brumfiel, A. Medina-Mardones, and J. Morgan. "A cochain
        level proof of Adem relations in the mod 2 Steenrod algebra."

        """

        def check_input(q, p):
            if not all(isinstance(i, int) for i in (q, p)) or \
                    any(isinstance(i, bool) for i in (q, p)):
                raise TypeError('q and p must be integers')
            if p < 0 or q < 0:
                raise ValueError('q and p must be non-negative')
            if q <= p:
                raise ValueError('q must be greater than p')

        def x(degree, action=0):
            transposition = (2, 1)
            key = tuple((1, 2) if (action + i) % 2 == 0 else transposition
                        for i in range(degree + 1))
            return BarrattEcclesElement({key: 1}, torsion=2)

        def action(simplex):
            if len(simplex) == 0:
                raise ValueError('simplex must be non-empty')
            identity, transposition = (1, 2), (2, 1)
            if simplex[0] == identity:
                start = 0
            elif simplex[0] == transposition:
                start = 1
            else:
                raise ValueError('defined for arity-2 simplices')
            for i, permutation in enumerate(simplex):
                expected = identity if (start + i) % 2 == 0 else transposition
                if permutation != expected:
                    raise ValueError('simplex is not alternating')
            return start

        def t_aw(left, right):
            answer = FreeModuleElement(torsion=2)
            for (left_key, left_coeff), (right_key, right_coeff) in product(
                    left.items(), right.items()):
                left_action = action(left_key)
                right_action = action(right_key)
                left_degree = len(left_key) - 1
                right_degree = len(right_key) - 1
                total_degree = left_degree + right_degree
                for epsilon in (0, 1):
                    for middle_degree in range(total_degree + 1):
                        final_degree = total_degree - middle_degree
                        if epsilon == 0:
                            coeff = comb(middle_degree + 1, left_degree + 1)
                        else:
                            coeff = comb(middle_degree, left_degree + 1)
                        coeff *= comb(final_degree, left_degree)
                        if coeff % 2 == 0:
                            continue
                        first = x(middle_degree, right_action)
                        second = x(final_degree, right_action + epsilon)
                        if left_action % 2:
                            first, second = second, first
                        answer += answer.create({
                            (next(iter(first)), next(iter(second))):
                            left_coeff * right_coeff
                        })
            return answer

        def phi(left, right):
            answer = FreeModuleElement(torsion=2)
            for (front, back), diagonal_coeff in left.diagonal().items():
                front = BarrattEcclesElement({front: 1}, torsion=2)
                back = BarrattEcclesElement({back: 1}, torsion=2)
                for (middle, final), coeff in t_aw(back, right).items():
                    answer += answer.create({
                        (next(iter(front)), middle, final):
                        diagonal_coeff * coeff
                    })
            return answer

        def compose_triple(outer, left, right):
            return outer.compose(left, 1).compose(right, 3)

        def operadic_phi(left, right):
            answer = BarrattEcclesElement(torsion=2)
            for (outer, middle, final), coeff in phi(left, right).items():
                outer = BarrattEcclesElement({outer: 1}, torsion=2)
                middle = BarrattEcclesElement({middle: 1}, torsion=2)
                final = BarrattEcclesElement({final: 1}, torsion=2)
                answer += coeff * compose_triple(outer, middle, final)
            return answer

        def ez_paths(rows, cols):
            if (rows, cols) == (0, 0):
                return [((0, 0),)]
            answer = []
            if rows > 0:
                for path in ez_paths(rows - 1, cols):
                    answer.append(path + ((rows, cols),))
            if cols > 0:
                for path in ez_paths(rows, cols - 1):
                    answer.append(path + ((rows, cols),))
            return answer

        def pair_product(first_degree, second_degree):
            first = next(iter(x(first_degree)))
            second = next(iter(x(second_degree)))
            answer = []
            for path in ez_paths(first_degree, second_degree):
                answer.append(tuple((first[i], second[j]) for i, j in path))
            return answer

        def pair_to_sigma4(pair):
            outer_action = 0 if pair[0] == (1, 2) else 1
            middle_action = 0 if pair[1] == (1, 2) else 1
            # The inclusion V_4 -> D_8 used in the paper sends the second
            # generator to bc, hence the second and third D_8 coordinates
            # agree.
            outer = x(0, outer_action)
            middle = x(0, middle_action)
            final = x(0, middle_action)
            return next(iter(compose_triple(outer, middle, final)))[0]

        def aw_pair_simplex(simplex):
            answer = FreeModuleElement(torsion=2)
            for i in range(len(simplex)):
                left_key = tuple(pair[0] for pair in simplex[:i + 1])
                right_key = tuple(pair[1] for pair in simplex[i:])
                left = BarrattEcclesElement({left_key: 1}, torsion=2)
                right = BarrattEcclesElement({right_key: 1}, torsion=2)
                left.preferred_rep()
                right.preferred_rep()
                if left and right:
                    answer += answer.create({
                        (next(iter(left)), next(iter(right))): 1
                    })
            return answer

        def psi_pair_simplex(simplex):
            answer = BarrattEcclesElement(torsion=2)
            for (left, right), coeff in aw_pair_simplex(simplex).items():
                left = BarrattEcclesElement({left: 1}, torsion=2)
                right = BarrattEcclesElement({right: 1}, torsion=2)
                answer += coeff * operadic_phi(left, right)
            return answer

        def join(prefix, element):
            answer = BarrattEcclesElement(torsion=2)
            for key, coeff in element.items():
                answer += answer.create({prefix + key: coeff})
            answer.preferred_rep()
            return answer

        def j_psi(first_degree, second_degree):
            answer = BarrattEcclesElement(torsion=2)
            for simplex in pair_product(first_degree, second_degree):
                for i in range(len(simplex)):
                    prefix = tuple(pair_to_sigma4(pair)
                                   for pair in simplex[:i + 1])
                    answer += join(prefix, psi_pair_simplex(simplex[i:]))
            return answer

        def j_transposition_23(first_degree, second_degree):
            transposition = SymmetricGroupElement((1, 3, 2, 4))
            answer = BarrattEcclesElement(torsion=2)
            for simplex in pair_product(first_degree, second_degree):
                simplex = tuple(pair_to_sigma4(pair) for pair in simplex)
                for i in range(len(simplex)):
                    moved = tuple(tuple(SymmetricGroupElement(sigma) *
                                        transposition)
                                  for sigma in simplex[i:])
                    answer += answer.create({simplex[:i + 1] + moved: 1})
            answer.preferred_rep()
            return answer

        check_input(q, p)
        answer = j_psi(q, p) + j_psi(p, q) + j_transposition_23(q, p)
        answer.preferred_rep()
        return answer

    @staticmethod
    def adem_coboundary(q, p):
        r"""Complete cochain-level Adem coboundary at the even prime.

        For integers :math:`q > p \geq 0`, this method returns

        .. math::

           A_{q,p} =
           \widetilde A_{q,p} + \widehat{NS}_{q,p}
           + \widehat{NS}_{p,q} \in \mathcal E(4),

        where :math:`\widetilde A_{q,p}` is
        :meth:`adem_homotopy` and

        .. math::

           \widehat{NS}_{q,p} =
           \sum_{\substack{0<a\leq\ell\\a\equiv\ell\pmod{\mathbb Z}}}
           \binom{p-\ell-a}{p-2\ell}
           \binom{p-\ell+a}{p-2\ell}
           \widetilde x_{q-p+2\ell+1}
           \circ
           (\widetilde x_{p-\ell-a},
            \widetilde x_{p-\ell+a}).

        All coefficients are reduced modulo 2. The action of this complete
        element on :math:`\alpha^{\otimes 4}`, rather than the action of
        :meth:`adem_homotopy` alone, is the cochain whose coboundary is the
        symmetric Adem sum in [BMM, Theorem 0.2.5].

        PARAMETERS
        ----------
        q, p : :class:`int`
            Non-negative integers with ``q > p``.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The arity-4, degree-:math:`p+q+1` complete Adem coboundary
            element over :math:`\mathbb F_2`.

        REFERENCES
        ----------
        [BMM]: G. Brumfiel, A. Medina-Mardones, and J. Morgan. "A cochain
        level proof of Adem relations in the mod 2 Steenrod algebra."

        """

        def non_symmetric_correction(outer_index, inner_index):
            answer = BarrattEcclesElement(torsion=2)

            # twice_l and twice_a encode the half-integers ell and a.
            # The binomial lower index forces twice_l <= inner_index.
            for twice_l in range(1, inner_index + 1):
                first_twice_a = 1 if twice_l % 2 else 2
                for twice_a in range(first_twice_a, twice_l + 1, 2):
                    lower = inner_index - twice_l
                    left_degree = (
                        2 * inner_index - twice_l - twice_a) // 2
                    right_degree = (
                        2 * inner_index - twice_l + twice_a) // 2
                    coefficient = (
                        comb(left_degree, lower)
                        * comb(right_degree, lower)
                    ) % 2
                    if coefficient == 0:
                        continue

                    outer_degree = (
                        outer_index - inner_index + twice_l + 1)
                    if outer_degree < 0:
                        continue

                    outer = BarrattEccles.may_steenrod_map(
                        2, outer_degree, torsion=2)
                    left = BarrattEccles.may_steenrod_map(
                        2, left_degree, torsion=2)
                    right = BarrattEccles.may_steenrod_map(
                        2, right_degree, torsion=2)
                    answer += outer.compose(left, 1).compose(right, 3)

            return answer

        answer = BarrattEccles.adem_homotopy(q, p)
        answer += non_symmetric_correction(q, p)
        answer += non_symmetric_correction(p, q)
        answer.preferred_rep()
        return answer
