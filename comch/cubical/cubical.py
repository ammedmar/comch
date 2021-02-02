from ..free_module import FreeModuleElement
from ..symmetric import SymmetricRingElement
from ..utils import pairwise
from itertools import combinations, product


class Cube(tuple):
    r"""A cube :math:`(I_1, \dots, I_n)`.

    A cube is a finite tuple of elements in :math:`\{0,1,2\}`, where we think
    of :math:`2` as the interval :math:`[0,1]` and :math:`0,1` as its endpoints.
    We identity these with faces of the infinite cube :math:`\mathbb I^\infty`.

    """

    def __init__(self, iterable):
        """Initializes *self*.

        PARAMETERS
        ----------
        interable : :class:'iterable'
            Used to create a :class:`tuple` of :class:`int` with values
            0, 1 or 2.

        EXAMPLE
        -------
        >>> print(Cube((1,2,0,2)))
        (1,2,0,2)

        """
        tuple.__init__(iterable)

    def __str__(self):
        return super.__str__(self).replace(', ', ',')

    @property
    def dimension(self):
        """The dimension of *self*.

        Defined as the number of values in the tuple that are equal to 2.

        RETURNS
        -------
        :class:`int`
            The dimension of *self*.

        EXAMPLE
        -------
        >>> Cube((1,2,0,2)).dimension
        2

        """
        return self.count(2)

    @property
    def intervals(self):
        """The positions of intervals in *self*.

        Corresponds to the tuple of indices where *self* contains 2.

        RETURNS
        -------
        :class:`tuple`
            The indices of intervals in *self*.

        EXAMPLE
        -------
        >>> Cube((1,2,0,2)).intervals
        (1, 3)

        """
        return tuple(idx for idx, x in enumerate(self) if x == 2)

    def face(self, i, epsilon):
        r"""The i-th :math:`\epsilon` face of *self*.

        Obtained by replacing the i-th entry of *self* by :math:`\epsilon`.

        RETURNS
        -------
        :class:`comch.simplicial.Simplex`
            The i-th face of *self*.

        EXAMPLE
        -------
        >>> Cube((1,2,0,2)).face(1, 0)
        (1, 2, 0, 1)

        """
        idx = self.intervals[i]
        answer = self[:idx] + ((epsilon + 1) % 2,) + self[idx + 1:]
        return Cube(answer)


class CubicalElement(FreeModuleElement):
    r"""Elements in an iterated tensor product of the chains on the
    infinite cube.

    The chains on the infinite cube :math:`C = C_\bullet(\Delta^\infty; R)`
    is the differential graded module :math:`C` with degree-:math:`n` part
    :math:`C_n` freely generated as an :math:`R`-module by cubes of
    dimension :math:`n`, and differential on these given by the sum of its
    faces with alternating signs. Explicitly, for $x \in C_n$ we have

    .. math::
        \partial x =
        \sum_{i = 1}^{n} (-1)^{i-1}(d^0_i x - d^1_i x)

    The differential graded module :math:`C` is isomorphic to

    .. math::
        \bigoplus_{k \geq 0} C_\bullet^{CW}(I; R)^{\otimes k}

    where :math:`I` is the interval with its usual cellular structure.

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

    """

    def __init__(self, data=None, torsion=None):
        """Initializes *self*.

        PARAMETERS
        ----------
        data : :class:`int` or ``None``, default: ``None``
            Dictionary representing a linear cobination of basis elements.
            Items in the dictionary correspond to `basis element: coefficient`
            pairs. Each basis element must create a :class:`tuple` of
            :class:`comch.cubical.Cube` and `coefficient` must be an
            :class:`int`.
        torsion : :class:`int` positive or :class:`string` equal to 'free'.
            The torsion of the underlying ring.

        EXAMPLE
        -------
        >>> x = CubicalElement({((0,2), (0,1)): 1,\
                                ((2,1), (1,2)): -1,\
                                ((1,2), (2,0)): 1})
        >>> print(x)
        ((0,2),(0,1)) - ((2,1),(1,2)) + ((1,2),(2,0))

        """

        if data:
            new_data = {}
            for k, v in data.items():
                new_k = tuple(Cube(cube) for cube in k)
                new_data[new_k] = v
            data = new_data

        super(CubicalElement, self).__init__(
            data=data, torsion=torsion)

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

        >>> x = CubicalElement({((2,), (1,)): 1, ((0,), (2,)): 1})
        >>> print(x._latex_())
        [01] \otimes [1] + [0] \otimes [01]

        """
        string = str(self)
        string = string.replace('1', '[1]').replace('0', '[0]')
        string = string.replace('2,', '[01],').replace(',2', ',[01]')
        string = string.replace(',),(', r' \otimes ')
        string = string.replace('),(', r' \otimes ')
        string = string.replace(')', '').replace('((', '')
        string = string.replace(')', '').replace('))', '')
        string = string.replace(',', '')

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
        >>> x = CubicalElement({((0,2), (0,1)): 1})
        >>> x.arity
        2

        """
        arities = set(len(multicube) for multicube in self.keys())
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
        >>> x = CubicalElement({((0,2), (0,1)): 1})
        >>> x.degree
        1

        """
        degrees = {sum(cube.dimension for cube in k) for k in self.keys()}
        if len(degrees) != 1:
            return None
        return degrees.pop()

    def boundary(self):
        """Boundary of *self*.

        As defined in the class's docstring.

        RETURNS
        _______
        :class:`comch.cubical.CubicalElement`
            The boundary of *self* as an element in a tensor product of
            differential graded modules.

        EXAMPLE
        -------
        >>> x = CubicalElement({((0, 2), (2, 1)): 1})
        >>> print(x.boundary())
        ((0,1),(2,1)) - ((0,0),(2,1)) - ((0,2),(1,1)) + ((0,2),(0,1))

        """
        answer = self.zero()
        for k, v in self.items():
            for idx, cube in enumerate(k):
                acc_dim = sum((cube_l.dimension for cube_l in k[:idx]))
                for i in range(cube.dimension):
                    for epsilon in (0, 1):
                        new_cube = cube.face(i, epsilon)
                        new_k = k[:idx] + (new_cube,) + k[idx + 1:]
                        sign_exp = (acc_dim + i + epsilon) % 2
                        answer += answer.create({new_k: v * (-1) ** sign_exp})
        return answer

    def __rmul__(self, other):
        """Left action: *other* ``*`` *self*

        Left multiplication by a symmetric group element or an integer.
        Defined up to signs on basis elements by permuting the tensor factor.

        PARAMETERS
        ----------
        other : :class:`int` or :class:`comch.symmetric.SymmetricElement`.
            The symmetric ring element left acting on *self*.

        RETURNS
        _______
        :class:`comch.cubical.CubicalElement`
            The product: *other* ``*`` *self* with Koszul's sign convention.

        EXAMPLE
        -------
        >>> x = CubicalElement({((0, 2), (1, 2)): 1})
        >>> t = SymmetricRingElement({(2, 1): 1})
        >>> print(t * x)
        - ((1,2),(0,2))
        >>> print(3 * x)
        3((0,2),(1,2))

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

        def sign(perm, multicube):
            weights = [cube.dimension % 2 for cube in multicube]
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
        r"""Iterated Serre diagonal applied at a specific tensor factor.

        The Serre diagonal is the chain map :math:`\Delta \colon C \to C \otimes C`
        defined on the chains of the infinite cube by the formula

        .. math::
            \Delta (x_1 \otimes \cdots \otimes x_n) =
            \sum \pm \left( x_1^{(1)} \otimes \cdots \otimes x_n^{(1)} \right) \otimes
            \left( x_1^{(2)} \otimes \cdots \otimes x_n^{(2)} \right),

        where the sign is determined using the Koszul convention, and we are
        using Sweedler’s notation

        .. math::
            \Delta(x_i) = \sum x_i^{(1)} \otimes x_i^{(2)}.

        It is coassociative, :math:`(\Delta \otimes \mathrm{id}) \Delta =
        (\mathrm{id} \otimes \Delta) \Delta`, so it has a well defined iteration
        :math:`\Delta^k`, and for every :math:`i \in \{1, \dots, r\}`, there is map
        :math:`C^{\otimes r} \to C^{\otimes k+r}` sending
        :math:`(x_1 \otimes \cdots \otimes x_n)` to
        :math:`(x_1 \otimes \cdots \otimes \Delta^k(k_i) \cdots \otimes x_n)`.

        PARAMETERS
        ----------
        times : :class:`int`
            The number of times the Serre diagonal is composed with itself.
        coord : :class:`int`
            The tensor position on which the iterated diagonal acts.

        RETURNS
        _______
        :class:`comch.cubical.CubicalElement`
            The action of the iterated Serre diagonal on *self*.

        EXAMPLE
        -------
        >>> x = CubicalElement({((2,2), (2,)):1})
        >>> print(x.iterated_diagonal(1,2))
        ((2,2),(2,),(1,)) + ((2,2),(0,),(2,))

        """

        def sign(p):
            """Counts the number of pairs appearing in reversed order."""
            to_count = filter(lambda x: x[0] > x[1], combinations(p, 2))
            sign_exp = sum(1 for _ in to_count) % 2
            return (-1) ** sign_exp

        def elementary_summand(fixed, i):
            """Models as a function the element 0,...,0,2,1,...,1 appearing
            as one of the summands of the iterated diagonal of an interval."""
            if i < fixed:
                return 0
            elif i == fixed:
                return 2
            else:
                return 1

        if self.degree is None:
            raise TypeError(f'only for homogeneous elements')
        if self.arity < coord:
            raise TypeError(f'arity = {self.arity} < coord = {coord}')

        answer = self.zero()
        for k, v in self.items():
            left, cube, right = k[:coord - 1], k[coord - 1], k[coord:]
            intervals = cube.intervals
            base = [i for idx, i in enumerate(cube) if idx not in intervals]
            for p in product(range(times + 1), repeat=cube.count(2)):
                multibase = [list(base) for _ in range(times + 1)]
                for idx, fixed in enumerate(p):
                    at = intervals[idx]
                    for i, new_base in enumerate(multibase):
                        to_insert = elementary_summand(fixed, i)
                        new_base.insert(at, to_insert)
                new_k = tuple(Cube(x) for x in multibase)
                answer += answer.create({left + new_k + right: v * sign(p)})
        return answer

    def join(self):
        r"""Join of *self*.

        The join is a map from :math:`C_\bullet^{\otimes 2}`
        to :math:`C_\bullet` of degree 1. It is define by

        .. math::
           \begin{aligned}
           (x_1 \otimes \cdots \otimes x_n) \ast (y_1 \otimes \cdots \otimes y_n)
           = (-1)^{|x|} \sum_{i=1}^n x_{<i} \epsilon(y_{<i}) \otimes
           x_i \ast y_i \otimes \epsilon(x_{>i})y_{>i},
           \end{aligned}

        where

        .. math::
           \begin{aligned}
           x_{<i} & = x_1 \otimes \cdots \otimes x_{i-1}, &
           y_{<i} & = y_1 \otimes \cdots \otimes y_{i-1}, \\
           x_{>i} & = x_{i+1} \otimes \cdots \otimes x_n, &
           y_{>i} & = y_{i+1} \otimes \cdots \otimes y_n,
           \end{aligned}

        with the convention

        .. math:: x_{<1} = y_{<1} = x_{>n} = y_{>n} = 1 \in \mathbb Z,

        and the only non-zero values of :math:`x_i \ast y_i` are

        .. math:: \ast([0] \otimes [1]) = [0, 1], \qquad  \ast([1] \otimes [0]) = -[0, 1].

        PARAMETERS
        ----------
        :class:`comch.cubical.CubicalElement`, of arity 2
            The element to take the join of.

        RETURNS
        _______
        :class:`comch.cubical.CubicalElement`
            The join of *self*.

        Examples
        --------
        >>> x = CubicalElement({((0, 0, 1), (1, 0, 0)): 1})
        >>> print(x.join())
        ((2,0,0),) - ((0,0,2),)

        """

        def is_zero(left, right):
            """Two conditions need to be satisfied for a triple
            (cube1, cube2, i) give a nonzero i-join: no intervals
            in cube1 can have indices greater or equal to i, and
            no intervals in cube2 can have indices less than or
            equal to i."""
            if left == tuple() or right == tuple():
                return False
            if isinstance(right, int):
                return right <= max(left)
            if isinstance(left, int):
                return left >= min(right)

        def _join(i, cube1, cube2, sign_exp):
            """the i-th elementary join keeping track of signs."""
            sign_exp += cube1.dimension
            cube = Cube(cube1[:i] + (2,) + cube2[i + 1:])
            p, q = cube1[i], cube2[i]
            if (p, q) == (0, 1):
                return cube, sign_exp % 2
            elif (p, q) == (1, 0):
                return cube, (1 + sign_exp) % 2
            else:
                return None, None

        if not self:
            return self

        if self.degree is None:
            raise TypeError(f'only for homogeneous elements')

        answer = self.zero()
        for k, v in self.items():
            for indices in combinations(range(len(k[0])), self.arity - 1):
                skip = False
                for i, (cube1, cube2) in zip(indices, pairwise(k)):
                    if (is_zero(cube1.intervals, i) or
                            is_zero(i, cube2.intervals)):
                        skip = True
                        break
                if not skip:
                    non_zero = True
                    sign_exp = 0
                    cube = k[0]
                    for i, next_cube in zip(indices, k[1:]):
                        cube, sign_exp = _join(i, cube, next_cube, sign_exp)
                        if cube is None:
                            non_zero = False
                            break
                    if non_zero:
                        answer += answer.create({(cube,): (-1) ** sign_exp * v})
        return answer


class Cubical:
    """Produces cubical elements of interest."""

    @staticmethod
    def standard_element(n, times=1, torsion=None):
        r"""The chain represented by the cube :math:`[0,1]^{n}`.

        PARAMETERS
        ----------
        n : :class:`int`
            The dimension of the standard cube considered.
        times : :class:`int`
            The number of tensor copies.
        torsion : :class:`int` positive or :class:`string` equal to 'free'
        The torsion of the underlying ring.

        EXAMPLES
        --------
        >>> print(Cubical.standard_element(3, 2))
        ((2,2,2),(2,2,2))
        """
        key = ((2,) * n,) * times
        return CubicalElement({key: 1}, torsion=torsion)
