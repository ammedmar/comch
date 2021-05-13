from ..free_module import FreeModuleElement

from ..symmetric import SymmetricGroupElement
from ..symmetric import SymmetricRingElement, SymmetricRing

from ..surjection import SurjectionElement
from ..utils import partitions, pairwise
from itertools import product, combinations
from math import floor, factorial


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
        """Initializes *self*.

        PARAMETERS
        ----------
        data : :class:`dict` or ``None``, default: ``None``
            Dictionary representing a linear combination of basis elements. Items
            in the dictionary correspond with pairs `basis_element: coefficient`.
            Each basis_element must create a :class:`tuple` of
            :class:`commch.symmetric.SymmetricGroupElement` and `coefficient` must
            be an :class:`int`.
        torsion : :class:`int` or :class:`string` 'free', default 'free'
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
        the complexity of an element agrees with its degree. It is proven in
        [BF] that the subcomplex generated by basis Barratt-Eccles elements of
        complexity at most :math:`n` define a suboperad of :math:`\mathcal E`
        modeling an :math:`E_{n+1}`-operad.

        RETURNS
        _______
        :class:`int`
            The complexity of *self*.

        EXAMPLE
        -------
        >>> BarrattEcclesElement({((1,2,3), (1,3,2), (1,2,3)): 1}).complexity
        1

        """
        complexity = 0
        for k in self.keys():
            for i, j in combinations(range(1, max(k[0]) + 1), 2):
                seq = (filter(lambda x: x == i or x == j, perm) for perm in k)
                cpxty = len([p for p, q in pairwise(seq) if p != q]) - 1
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

        This operad map preserves the :math:`E_n`-filtration.

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
    def may_steenrod_structure(arity, degree, torsion=None):
        r"""Representative of the requested Steenrod product.

        Returns the image under the Steenrod-Adem structure constructed in
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
            The torsion of the underlying ring.

        RETURNS
        _______
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The image of the basis element under the Steenrod-Adem structure.

        REFERENCES
        ----------
        [KMM]: Kaufmann, R. M., & Medina-Mardones, A. M. (2020). Chain level
        Steenrod operations. arXiv preprint arXiv:2010.02571.
        
        SEE ALSO
        --------
        :meth:`comch.surjection.Surjection.may_steenrod_structure`.

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

    @staticmethod
    def steenrod_operation(p, s, q, bockstein=False):
        r"""Chain level representative of :math:`P_s` or :math:`\beta P_s`.

        Where

        .. math::
            P_s : H_\bullet(A; \mathbb F_2)
            \to
            H_{\bullet + s}(A; \mathbb F_2)

        and, for :math:`p > 2`,

        .. math::
           P_s : H_\bullet(A; \mathbb F_p)
           \to
           H_{\bullet + 2s(p-1)}(A; \mathbb F_p)
           \qquad
           \beta P_s : H_\bullet(A; \mathbb F_p)
           \to
           H_{\bullet+2s(p-1)-1}(A; \mathbb F_p)

        are the Steenrod operations.

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

        SEE ALSO
        --------
        :meth:`comch.surjection.Surjection.steenrod_operation`

        """

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
                return BarrattEcclesElement(torsion=2)
        else:
            b = int(bockstein)
            # Serre convention: v(2j)=(-1)^j & v(2j+1)=v(2j)*m! w/ m=(p-1)/2
            coeff = (-1) ** (floor(q / 2) + s)
            if q / 2 - floor(q / 2):
                coeff *= factorial((p - 1) / 2)
            # degree of the element: (2s-q)(p-1)-b
            d = (2 * s - q) * (p - 1) - b
            if d < 0:
                return BarrattEcclesElement(torsion=p)

        return int(coeff) * BarrattEccles.may_steenrod_structure(
            p, d, torsion=p)
