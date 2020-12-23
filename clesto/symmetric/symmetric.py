from clesto.module import Module_element
from itertools import product


class SymmetricGroup_element(tuple):
    r"""Element in a finite symmetric group.

    We refer to elements in the group of permutations of :math:`r` elements
    :math:`\mathrm S_r` as symmetric elements of arity :math:`r`. An
    element :math:`\pi \in \mathrm S_r` can be thought of as a bijections
    from :math:`\{1,\dots,r\}` to itself, and can be represented by the
    tuple of its images :math:`(\pi(1), \dots, \pi(r))`.

    PARAMETERS
    ----------
        iterable : 'iterable'
            Iterable representing a permutation of `(1,...,r)`.

    EXAMPLE
    -------

    >>> print(SymmetricGroup_element((1,3,2)))
    (1,3,2)

    """

    def __str__(self):
        s = super().__str__()
        return s.replace(', ', ',')

    @property
    def sign(self):
        """Returns the sign of *self*.

        The sign is defined as the mod 2 number of transpositions required
        to express the element.

        RETURNS
        -------

        :class: `int`
            The sign of *self*.

        EXAMPLE
        -------

        >>> SymmetricGroup_element((5,2,4,3,1)).sign
        1

        """
        if set(self) != set(range(1, len(self) + 1)):
            raise TypeError(f'defined for permutations of (1,...,r) ' +
                            f'only not {self}')
        cycles = self.to_cycles()
        return (-1)**sum(len(cycle) - 1 for cycle in cycles)

    def to_cycles(self, singletons=False):
        """Transforms from bijection to collection of cycles.

        RETURNS
        -------

        :class: `list` of `tuple`
            The representation of *self* as a product of cycles.

        EXAMPLE
        -------

        >>> SymmetricGroup_element((5,2,4,3,1)).to_cycles()
        [(1, 5), (3, 4)]

        """
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

    @property
    def arity(self):
        """Arity of *self*

        The arity of a symmetric group element is defined as the
        cardinality of its domain.

        RETURNS
        -------

        :class: `list` of 'tuple'
            The arity of *self*.

        EXAMPLE
        -------

        >>> SymmetricGroup_element((5,2,4,3,1)).arity
        5

        """
        return max(self)

    def __mul__(self, other):
        r"""Product: *self* * *other*.

        This product agrees with the composition of bijections:
        *self* :math:`\circ` *other*.

        RETURNS
        -------

        :class: `clesto.symmetric.symmetric.SymmetricGroup_element` object
            The product of *self* and *other*.

        EXAMPLE
        -------

        >>> x = SymmetricGroup_element((1,3,2))
        >>> y = SymmetricGroup_element((2,3,1))
        >>> print(y * x)
        (2,1,3)

        """
        # default to method in *other*
        if not isinstance(other, SymmetricGroup_element):
            self = SymmetricRing_element({self: 1}, torsion=other.torsion)
            return other.__rmul__(self)

        if self.arity != other.arity:
            raise TypeError('Unequal arity attribute')

        return SymmetricGroup_element(tuple(self[i - 1] for i in other))

    def inverse(self):
        """Multiplicative inverse: *self*:math:`^{-1}`.

        RETURNS
        -------

        :class: `clesto.symmetric.symmetric.SymmetricGroup_element` object
            The multiplicative inverse of *self*.

        EXAMPLES
        -------

        >>> pi = SymmetricGroup_element((2,3,1))
        >>> print(pi.inverse())
        (3,1,2)

        """
        inverse = tuple(self.index(i + 1) + 1 for i in range(self.arity))
        return SymmetricGroup_element(inverse)

    def __pow__(self, times):
        """Iterated product of *self*: *self* * ... * *self*.

        RETURNS
        -------

        :class: `clesto.symmetric.symmetric.SymmetricGroup_element` object
            The product of *self* with itself *times* number of times.

        EXAMPLE
        -------

        >>> x = SymmetricGroup_element((2,3,4,5,1))
        >>> print(x**5)
        (1,2,3,4,5)

        """
        if times == 0:
            return SymmetricGroup_element(tuple(range(1, max(self) + 1)))
        answer = self
        for i in range(times - 1):
            answer *= self

        return answer

    def compose(self, other, position):
        r"""Operadic compositions: *self* :math:`\circ_{position}` *other*.

        The (operadic) composition of symmetric elements is defined as follows:
        Given :math:`\pi \in \Sigma_r`, :math:`\tau \in \Sigma_{s}` and
        :math:`i \in \{1, \dots, r\}` produces a permutation in
        :math:`\Sigma_{r + s - 1}`. We begin by considering
        :math:`\{1, 2, \dots, r + s - 1\}` as an ordered set :math:`R` with
        :math:`r` elements by grouping the subset
        :math:`S = \{i, \dots, i+s-1\}` into a single element, then applying
        :math:`\pi` to :math:`R` and :math:`\sigma` to the :math:`S`, and,
        finally, forgetting the grupping. More precisely, for integers
        :math:`r, s \geq 1` and :math:`i \in \{1, \ldots, r\}` the partial
        composition is the linear map

        .. math:: \circ_i : \Sigma_r \otimes \Sigma_s \to \Sigma_{r+s-1}

        is defined for :math:`\pi = (\pi(1), \dots, \pi(r))` and
        :math:`\sigma = (\sigma(1), \dots, \sigma(s))` to be the sequence
        obtained by replacing in position :math:`i` of the sequence :math:`\pi`
        the sequence obtained by adding :math:`i-1` to the entries of :math:`s`
        and adding :math:`s-1` to the entries of :math:`\pi` that whose value
        is greater than :math:`i`.

        RETURNS
        -------

        :class: `clesto.symmetric.symmetric.SymmetricGroup_element` object
            The composition of *self* and *other*.

        EXAMPLE
        -------

        >>> x = SymmetricGroup_element((1,3,2))
        >>> y = SymmetricGroup_element((2,1))
        >>> print(x.compose(y, 1))
        (2,1,4,3)

        """
        s = len(other) - 1
        to_insert = tuple(i + position - 1 for i in other)
        at = self.index(position)
        shift = tuple(map(lambda i: i + s if i > position else i, self))
        answer = shift[:at] + to_insert + shift[at + 1:]
        return SymmetricGroup_element(answer)


class SymmetricRing_element(Module_element):
    """Elements in the (modular) integral group ring of finite symmetric groups.

    """

    def __init__(self, data=None, torsion=None):
        if data:
            data = {SymmetricGroup_element(k): v for k, v in data.items()}

        super(SymmetricRing_element, self).__init__(data=data,
                                                    torsion=torsion)

    def __str__(self):
        s = super().__str__()
        return s.replace(', ', ',')

    @property
    def arity(self):
        """Return the arity of *self* if homogeneous and ``None`` otherwise.

        *self* is said to be homogeneous if all basis elements belong to the
        same arity.

        >>> SymmetricRing_element({(5,2,4,3,1): 1}).arity
        5
        >>> SymmetricRing_element({(2,3,1): 1, (1,2): 1}).arity

        """
        if not self:
            return None

        arities = set(max(k) for k in self.keys())
        if len(arities) > 1:
            return None

        return arities.pop()

    def __mul__(self, other):
        """Linear product in the symmetric group ring: *self* * *other*.

        >>> p = SymmetricRing_element({(4,3,2,1): 1, (1,2,3,4): 2})
        >>> print(3 * p)
        3(4,3,2,1) + 6(1,2,3,4)
        >>> q = SymmetricRing_element({(4,1,2,3): 1})
        >>> print(p * q)
        (1,4,3,2) + 2(4,1,2,3)

        """

        if isinstance(other, int):
            return super().__rmul__(other)

        if not isinstance(other, SymmetricRing_element):
            return other.__rmul__(self)

        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')

        answer = self.zero()
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            answer[tuple(k1[i - 1] for i in k2)] += v1 * v2

        answer._reduce_rep()
        return answer

    def __pow__(self, times):
        """Iterated product of *self*: *self* * ... * *self*.

        >>> p = SymmetricRing_element({(4,3,2,1): 1, (1,2,3,4): 2})
        >>> p ** 2
        SymmetricRing_element({(1, 2, 3, 4): 3})

        """
        if times == 0:
            return SymmetricRing.identity_element(self.arity, self.torsion)
        answer = self.zero()
        for k, v in self.items():
            answer += self.create({k**times: v})
        return answer

    def compose(self, other, position):
        """Linear operadic compositions: *self* o_position *other*.

        >>> x = SymmetricRing_element({(2,3,1): 1, (1,2,3): -1})
        >>> y = SymmetricRing_element({(2,1): 1, (1,2): 1})
        >>> print(x.compose(y, 2))
        (3,2,4,1) + (2,3,4,1) - (1,3,2,4) - (1,2,3,4)

        """
        if not self:
            return self.zero()

        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')

        answer = self.zero()
        for (k1, v1), (k2, v2) in product(self.items(), other.items()):
            new_k = k1.compose(k2, position)
            new_v = v1 * v2
            to_add = self.create({new_k: new_v})
            answer += to_add
        return answer


class SymmetricRing():
    """Produce special elements in the group ring of finite symmetric groups.

    """
    @staticmethod
    def identity_element(arity, torsion=None):
        """Return (1,...,r)."""
        identity = tuple(range(1, arity + 1))
        return SymmetricRing_element({identity: 1}, torsion=torsion)

    @staticmethod
    def rotation_element(arity, torsion=None):
        """Return (2,3,...,r-1,1)."""
        rho = tuple(range(2, arity + 1)) + (1,)
        return SymmetricRing_element({rho: 1}, torsion=torsion)

    @staticmethod
    def transposition_element(arity, torsion=None):
        """Return (2,3,...,r-1,1) - (1,2,...,r)."""
        rho = tuple(range(2, arity + 1)) + (1,)
        identity = tuple(range(1, arity + 1))
        return SymmetricRing_element({rho: 1, identity: -1},
                                     torsion=torsion)

    @staticmethod
    def norm_element(arity, torsion=None):
        """Return (2,3,...,r-1,1) ** 0 + ... + (2,3,...,r-1,1) ** r."""
        rho = SymmetricRing.rotation_element(arity, torsion=torsion)
        answer = SymmetricRing_element(torsion=torsion)
        for i in range(arity):
            answer += rho**i
        return answer
