from ..free_module import FreeModuleElement
from ..symmetric import SymmetricRingElement
from ..utils import pairwise
from itertools import chain, product, combinations_with_replacement


class Simplex(tuple):
    r"""A simplex :math:`(v_0, \dots, v_n)`.

    A simplex is a finite non-decreasing tuple of non-negative integers.
    We identity these with faces of the infinite simplex :math:`\Delta^\infty`.

    """

    def __init__(self, iterable):
        """Initializes *self*.

        PARAMETERS
        ----------
        interable : :class:'iterable'
            Used to create a :class:`tuple` of :class:`int`.

        EXAMPLE
        -------
        >>> print(Simplex((1,2,4)))
        (1,2,4)

        """
        tuple.__init__(iterable)

    def __str__(self):
        return super.__str__(self).replace(', ', ',')

    @property
    def dimension(self):
        """The dimension of *self*.

        Defined as the length of the tuple minus one.

        RETURNS
        -------
        :class:`int`
            The dimension of *self*.

        EXAMPLE
        -------
        >>> Simplex((1,3,4,5)).dimension
        3

        """
        return len(self) - 1

    def face(self, i):
        """The i-th face of *self*.

        Obtained by removing the i-th entry of *self*.

        RETURNS
        -------
        :class:`comch.simplicial.Simplex`
            The i-th face of *self*.

        EXAMPLE
        -------
        >>> Simplex((1,3,4,5)).face(2)
        (1, 3, 5)

        """
        return Simplex(self[:i] + self[i + 1:])

    def degeneracy(self, i):
        """The i-th degeneracy of *self*.

        Obtained by repeating the i-th entry of the tuple.

        RETURNS
        -------
        :class:`comch.simplicial.Simplex`
            The i-th face of *self*.

        EXAMPLE
        -------
        >>> Simplex((1,3,4,5)).degeneracy(2)
        (1, 3, 4, 4, 5)

        """

        return Simplex(self[:i + 1] + self[i:])

    def coface(self, i):
        """The i-th coface of *self*.

        Obtained by adding 1 to each j-th entries with j
        greater or equal to i.

        RETURNS
        -------
        :class:`comch.simplicial.Simplex`
            The i-th coface of *self*.

        EXAMPLE
        -------
        >>> Simplex((1,3,4,5)).coface(2)
        (1, 4, 5, 6)

        """

        def d_i(i, j): return j + 1 if j >= i else j

        return tuple(d_i(i, j) for j in self)

    def codegeneracy(self, i):
        """The i-th codegeneracy of *self*.

        Obtained by subtracting 1 from each j-th entries with j
        greater than i.

        RETURNS
        -------
        :class:`comch.simplicial.Simplex`
            The i-th codegeneracy of *self*.

        EXAMPLE
        -------
        >>> Simplex((1,3,4,5)).codegeneracy(2)
        (1, 2, 3, 4)

        """

        def s_i(i, j): return j - 1 if j > i else j

        return tuple(s_i(i, j) for j in self)

    def is_degenerate(self):
        """Returns ``True`` if *self* is degenerate and ``False`` if not.

        A simplex is degenerate if it is empty or if contains equal consecutive
        values.

        RETURNS
        -------
        :class:`bool`
            ``True`` if *self* is degenerate and ``False`` if not.

        EXAMPLE
        -------
        >>> Simplex(()).is_degenerate()
        True
        >>> Simplex((1,1,2)).is_degenerate()
        True

        """
        consecutive_values = any([i == j for i, j in pairwise(self)])
        empty_simplex = (self.dimension == -1)
        return empty_simplex or consecutive_values

    def is_nondegenerate(self):
        """Returns ``True`` if *self* is nondegenerate and ``False`` if not.

        A simplex is nondegenerate if it is not empty and contains no equal
        consecutive values.

        RETURNS
        -------
        :class:`bool`
            ``True`` if *self* is nondegenerate and ``False`` if not.

        EXAMPLE
        -------
        >>> Simplex((1,2,5)).is_nondegenerate()
        True

        """
        return not self.is_degenerate()


class SimplicialElement(FreeModuleElement):
    r"""Elements in an iterated tensor product of the chains on the
    infinite simplex.

    The chains on the infinite simplex :math:`C = C_\bullet(\Delta^\infty; R)`
    is the differential graded module :math:`C` with degree-:math:`n` part
    :math:`C_n` freely generated as an :math:`R`-module by simplices of
    dimension :math:`n`, and differential on these given by the sum of its
    faces with alternating signs. Explicitly,

    .. math::
        \partial (v_0, \dots, v_n) =
        \sum_{i=0}^n (v_0, \dots, \widehat{v}_i, \dots, v_d).

    The degree-:math:`n` part of the tensor product :math:`C^{\otimes r}`
    and its differential are recursively defined by

    .. math::
        (C^{\otimes r})_n = \bigoplus_{i+j=n} C_i \otimes (C^{\otimes r})_n

    and

    .. math::
        \partial(c_1 \otimes c_2 \otimes \cdots \otimes c_r) =
        (\partial c_1) \otimes c_2 \otimes \cdots \otimes c_r +
        (-1)^{|c_1|} c_1 \otimes \partial (c_2 \otimes \cdots \otimes c_r).

    ATTRIBUTES
    ----------
    torsion : :class:`int` positive or :class:`string` equal to 'free'.
        The torsion of the underlying ring.
    dimension : :class:`int` non-negative.
        NOT SURE IF NEEDED.

    """

    def __init__(self, data=None, dimension=None, torsion=None):
        """Initializes *self*.

        PARAMETERS
        ----------
        data : :class:`int` or ``None``, default: ``None``
            Dictionary representing a linear cobination of basis elements.
            Items in the dictionary correspond to `basis_element: coefficient`
            pairs. Each basis_element must create a :class:`tuple` of
            :class:`comch.simplicial.Simplex` and `coefficient` must be an
            :class:`int`.
        dimension : :class:`int`
            NOT SURE IF NEEDED.
        torsion : :class:`int` positive or :class:`string` equal to 'free'.
            The torsion of the underlying ring.

        EXAMPLE
        -------
        >>> x = SimplicialElement({((0,), (0, 1, 2)): 1,\
                                   ((0, 1), (1, 2)): -1,\
                                   ((0, 1, 2), (2,)): 1})
        >>> print(x)
        ((0,),(0,1,2)) - ((0,1),(1,2)) + ((0,1,2),(2,))

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

        super(SimplicialElement, self).__init__(data=data,
                                                torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def _latex_(self):
        r"""Representation in LaTex.

        RETURNS
        -------
        :class:`string`
            A LaTex friendly representation.

        EXAMPLE
        -------
        >>> x = SimplicialElement({((0,), (0, 1, 2)): 1})
        >>> print(x._latex_())
        [0] \otimes [0,1,2]

        """
        string = str(self)
        string = string.replace(',),(', r'] \otimes [')
        string = string.replace('),(', r'] \otimes [')
        string = string.replace('((', '[')
        string = string.replace(',))', ']').replace('),)', ']')
        string = string.replace('))', ']')

        return string

    @property
    def arity(self):
        """Arity of *self*.

        Defined as ``None`` if *self* is not homogeneous. The arity of a basis
        element is defined as the number of tensor factors making it.

        RETURNS
        -------
        :class:`int` positive or ``None``.
            The length of the keys of *self* or ``None`` if not well defined.

        EXAMPLE
        -------
        >>> x = SimplicialElement({((0,), (0, 1, 2)): 1})
        >>> x.arity
        2

        """
        arities = set(len(multispx) for multispx in self.keys())
        if len(arities) != 1:
            return None
        return arities.pop()

    @property
    def degree(self):
        """Degree of *self*.

        Defined as ``None`` if self is not homogeneous. The degree of a basis
        element agrees with the sum of the dimension of the simplices making it.

        RETURNS
        -------
        :class:`int` positive or ``None``.
            The sum of the dimensions of the simplices of every key of *self* or
            ``None`` if not well defined.

        EXAMPLE
        -------
        >>> x = SimplicialElement({((0,), (0, 1, 2)): 1})
        >>> x.degree
        2

        """
        degs = {sum(spx.dimension for spx in k) for k in self.keys()}
        if len(degs) != 1:
            return None
        return degs.pop()

    def boundary(self):
        """Boundary of *self*.

        As defined in the class's docstring.

        RETURNS
        _______
        :class:`comch.simplicical.SimplicialElement`
            The boundary of *self* as an element in a tensor product of
            differential graded modules.

        EXAMPLE
        -------
        >>> x = SimplicialElement({((0, 1), (1, 2)): 1})
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
                    answer += answer.create({new_k: v * (-1) ** sign_exp})
        return answer

    def __rmul__(self, other):
        """Left action: *other* ``*`` *self*

        Left multiplication by a symmetric group element or an integer.
        Defined up to signs on basis elements by permuting the tensor factor.

        PARAMETERS
        ----------
        other : :class:`int` or :class:`comch.simplicial.SimplicialElement`.
            The symmetric ring element left acting on *self*.

        RETURNS
        _______
        :class:`comch.simplicial.SimplicialElement`
            The product: *other* ``*`` *self* with Koszul's sign convention.

        EXAMPLE
        -------
        >>> x = SimplicialElement({((0, 1), (1, 2)): 1})
        >>> t = SymmetricRingElement({(2, 1): 1})
        >>> print(t * x)
        - ((1,2),(0,1))
        >>> print(3 * x)
        3((0,1),(1,2))

        """

        def check_input(self, other):
            """Symmetric ring element with same attributes."""
            if not isinstance(other, SymmetricRingElement):
                raise TypeError(f'__rmul__ by type int or \
                    SymmetricRingElement not {type(other)}')
            if self.torsion != other.torsion:
                raise TypeError('only defined for equal attribute torsion')
            if self.arity != other.arity:
                raise TypeError('Unequal arity attribute')

        def sign(perm, multispx):
            weights = [spx.dimension % 2 for spx in k1]
            sign_exp = 0
            for idx, i in enumerate(perm):
                right = [weights[perm.index(j)] for
                         j in perm[idx + 1:] if i > j]
                sign_exp += sum(right) * weights[idx]
            return (-1) ** (sign_exp % 2)

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

    def iterated_diagonal(self, times=1, coord=1):
        r"""Iterated Alexander-Whitney diagonal applied at a specific tensor factor.

        The AW diagonal is the chain map :math:`\Delta \colon C \to C \otimes C`
        defined on the chains of the infinite simplex by the formula

        .. math::
            \Delta((v_0, dots, v_n)) =
            \sum_{i=0}^n (v_0, dots, v_i) \otimes (v_i, dots, v_n).

        It is coassociative, :math:`(\Delta \otimes \mathrm{id}) \Delta =
        (\mathrm{id} \otimes \Delta) \Delta`, so it has a well defined iteration
        :math:`\Delta^k`, and for every :math:`i \in \{1, \dots, r\}`, there is map
        :math:`C^{\otimes r} \to C^{\otimes k+r}` sending
        :math:`(x_1 \otimes \cdots \otimes x_n)` to
        :math:`(x_1 \otimes \cdots \otimes \Delta^k(k_i) \cdots \otimes x_n)`.

        PARAMETERS
        ----------
        times : :class:`int`
            The number of times the AW diagonal is composed with itself.
        coord : :class:`int`
            The tensor position on which the iterated diagonal acts.

        RETURNS
        _______
        :class:`comch.simplicial.SimplicialElement`
            The action of the iterated AW diagonal on *self*.

        EXAMPLE
        -------
        >>> x = SimplicialElement({((0, 1, 2), ): 1})
        >>> print(x.iterated_diagonal())
        ((0,),(0,1,2)) + ((0,1),(1,2)) + ((0,1,2),(2,))

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

    def one_reduced(self):
        """Returns the 1-reduction of *self*.

        The 1-reduction map is the map induced by the collapse of the
        1-skeleton of the infinite simplex.

        RETURNS
        _______
        :class:`comch.simplicial.SimplicialElement`
            The preferred representative of *self*.

        EXAMPLE
        -------
        >>> x = SimplicialElement({((1,2), (2,3,4)): 1})
        >>> print(x.one_reduced())
        0

        """
        answer = self.zero()
        for k, v in self.items():
            if all(spx.dimension != 1 for spx in k):
                answer += self.create({k: v})
        return answer

    def preferred_rep(self):
        """Preferred representative of *self*.

        Removes pairs `basis element: coefficient` which satisfy either of:
        1) The basis element has a degenerate tensor factor, or 2) the
        coefficient is 0.

        RETURNS
        _______
        :class:`comch.simplicial.SimplicialElement`
            The preferred representative of *self*.

        EXAMPLE
        -------
        >>> print(SimplicialElement({((1,3), (1,1)): 1}))
        0

        """
        for k, v in self.items():
            if any([spx.is_degenerate() for spx in k]):
                self[k] = 0

        super().preferred_rep()


class Simplicial:
    """Produces simplicial elements of interest."""

    @staticmethod
    def standard_element(n, times=1, torsion=None):
        r"""The chain represented by the simplex :math:`(0, \dots, n)`.

        PARAMETERS
        ----------
        n : :class:`int`
            The dimension of the standard simplex considered.
        times : :class:`int`
            The number of tensor copies.
        torsion : :class:`int` positive or :class:`string` equal to 'free'
        The torsion of the underlying ring.

        EXAMPLES
        --------
        >>> print(Simplicial.standard_element(3, 2))
        ((0,1,2,3),(0,1,2,3))

        """
        key = (tuple(range(n + 1)),) * times
        return SimplicialElement({key: 1}, torsion=torsion)
