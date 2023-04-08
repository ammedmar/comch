from collections import Counter


class FreeModuleElement(Counter):
    r"""Element in a free :math:`\mathbb{Z}/n \mathbb{Z}`-module.

    Let :math:`R` be a ring and :math:`B` a set. The free :math:`R`-module
    generated by :math:`B` consists of all :math:`R`-linear combination of
    elements in :math:`B`

    .. math::
        R[B] = \Big\{ \sum_i r_ib_i\ |\ r_i \in R, b_i \in B \Big\}.

    This class models elements in free :math:`\mathbb Z/n \mathbb Z`-modules
    with :math:`\mathbb Z/0 \mathbb Z = \mathbb Z`. The ring is specified via
    the class attribute *torsion* corresponding to the non-negative integer
    :math:`n`.

    ATTRIBUTES
    ----------
    torsion : :class:`int`
        The non-neg int :math:`n` of the ring :math:`\mathbb Z/n \mathbb Z`.

    """

    default_torsion = 0
    """Class attribute: Used if :attr:`torsion` is ``None`` during
    initialization."""

    def __init__(self, data=None, torsion=None):
        """Initializes *self*.

        PARAMETERS
        ----------
        data : ``dict``
            Dictionary representing a linear combination of basis elements.
            Items in the dict correspond to pairs (basis_element: coefficient).
        torsion : :class:`int`
            The non-neg int :math:`n` of the ring :math:`\mathbb Z/n \mathbb Z`.

        EXAMPLE
        -------
        >>> print(FreeModuleElement())
        0
        >>> print(FreeModuleElement({'a': 1, 'b': -1, 'c': 0}))
        a - b

        """
        if torsion is None:
            torsion = type(self).default_torsion
        self.torsion = torsion
        super(FreeModuleElement, self).__init__(data)
        self.preferred_rep()

    def __hash__(self):
        return hash(frozenset(self))

    def __str__(self):
        """Coefficient first representation."""
        if not self:
            return '0'
        else:
            answer = ''
            for k, v in self.items():
                if v < -1:
                    answer += f'- {abs(v)}{str(k)} '
                elif v == -1:
                    answer += f'- {str(k)} '
                elif v == 1:
                    answer += f'+ {str(k)} '
                elif v > 1:
                    answer += f'+ {v}{str(k)} '
            if answer[0] == '+':
                answer = answer[2:]
            return answer[:-1]

    def __add__(self, other):
        """Addition: *self* + *other*.

        PARAMETERS
        ----------
        other : :class:`comch.free_module.FreeModuleElement` object
            The element to add to *self*.

        RETURNS
        -------
        :class:`comch.free_module.FreeModuleElement` object
            The sum of *self* and *other*.

        EXAMPLE
        -------
        >>> FreeModuleElement({'a': 1, 'b': 2}) + FreeModuleElement({'a': 1})
        FreeModuleElement({'a': 2, 'b': 2})

        """
        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')
        answer = self.create(self)
        answer.update(other)
        answer.preferred_rep()
        return answer

    def __sub__(self, other):
        """Diference: *self* - *other*.

        PARAMETERS
        ----------
        other : :class:`comch.free_module.FreeModuleElement` object
            The element to subtract from *self*.

        RETURNS
        -------
        :class:`comch.free_module.FreeModuleElement` object
            The difference of *self* and *other*.

        EXAMPLE
        -------
        >>> FreeModuleElement({'a': 1, 'b': 2}) - FreeModuleElement({'a': 1})
        FreeModuleElement({'b': 2})

        """
        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')
        answer = self.create(self)
        answer.subtract(other)
        answer.preferred_rep()
        return answer

    def __rmul__(self, c):
        """Left action: *c* * *self*.

        PARAMETERS
        ----------
        c : int
            The element to act *self* with.

        RETURNS
        -------
        :class:`comch.free_module.FreeModuleElement` object
            The action of *c* on *self*.

        EXAMPLE
        -------
        >>> 3 * FreeModuleElement({'a':1, 'b':2})
        FreeModuleElement({'b': 6, 'a': 3})

        """
        if not isinstance(c, int):
            raise TypeError(f'Act only by int not by type {type(c)}')

        scaled = {k: c * v for k, v in self.items()}
        return self.create(scaled)

    def __truediv__(self, c):
        """Division of *self* by *c*.

        When c has a multiplicative inverse in the ground ring, divides *self* by *c*.

        PARAMETERS
        ----------
        c : int
            The element to divide *self* by.

        RETURNS
        -------
        :class:`comch.free_module.FreeModuleElement` object
            The action of *1/c* on *self*.

        EXAMPLE
        -------
        >>> FreeModuleElement({'a': 1, 'b': 2}, torsion=5) / 3
        FreeModuleElement({'b': 4, 'a': 2})

        """
        if not isinstance(c, int):
            raise TypeError(f'Act only by int not by type {type(c)}')
        if self.torsion == 0:
            raise TypeError(f'Division defined for positive torsion only')

        inv = pow(c, -1, self.torsion)
        scaled = {k: inv * v for k, v in self.items()}
        return self.create(scaled)

    def __neg__(self):
        """Additive inverse: - *self*.

        RETURNS
        -------
        :class:`comch.free_module.FreeModuleElement` object
            The additive inverse of *self*.

        EXAMPLE
        -------
        >>> - FreeModuleElement({'a': 1, 'b': 2})
        FreeModuleElement({'a': -1, 'b': -2})

        """
        return self.__rmul__(-1)

    def __iadd__(self, other):
        """In place addition: *self* = *self* + *other*.

        PARAMETERS
        ----------
        other : :class:`comch.free_module.FreeModuleElement` object
            The element to add to *self*.

        EXAMPLE
        -------
        >>> x = FreeModuleElement({'a': 1, 'b': 2})
        >>> x += FreeModuleElement({'a': 3, 'b': 6})
        >>> x
        FreeModuleElement({'b': 8, 'a': 4})

        """
        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')
        self.update(other)
        self.preferred_rep()
        return self

    def __isub__(self, other):
        """In place difference: *self* = *self* - *other*.

        PARAMETERS
        ----------
        other : :class:`comch.free_module.FreeModuleElement` object
            The element to subtract from *self*.

        EXAMPLE
        -------
        >>> x = FreeModuleElement({'a': 1, 'b': 2})
        >>> x -= FreeModuleElement({'a': 3, 'b': 6})
        >>> x
        FreeModuleElement({'a': -2, 'b': -4})

        """
        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')
        self.subtract(other)
        self.preferred_rep()
        return self

    def preferred_rep(self):
        r"""The preferred representative of *self*.

        Consisting of pairs `basis_element: coefficient` with `coefficient`
        different from 0 and in the set :math:`\{1, \dots, r-1\}` if
        :attr:`torsion` is an :class:`int` denoted :math:`r`.

        EXAMPLE
        -------
        >>> FreeModuleElement({'a': 1, 'b': 2, 'c': 0})
        FreeModuleElement({'b': 2, 'a': 1})

        """
        # reducing coefficients mod torsion
        if self.torsion != 0:
            for k, v in self.items():
                self[k] = v % self.torsion

        # removing k:v pairs with v = 0
        zeros = [k for k, v in self.items() if not v]
        for k in zeros:
            del self[k]

    def set_torsion(self, torsion):
        """Sets the torsion of *self*.

        PARAMETERS
        ----------
        torsion : :class:`int`
            The new `torsion` of *self*

        EXAMPLE
        -------
        >>> FreeModuleElement({'a': 1, 'b': 2}).set_torsion(2)
        FreeModuleElement({'a': 1})

        """
        setattr(self, 'torsion', torsion)
        self.preferred_rep()
        return self

    def create(self, other=None):
        """Initializes with the same type and attribute values as *self*.

        PARAMETERS
        ----------
        other : dict or None, default: ``None``
            Data to be initialized.

        RETURNS
        -------
        type(*self*) object
            The initialized object with the given data

        EXAMPLE
        -------
        >>> x =  FreeModuleElement({'a': 1})
        >>> x + x.create({'b': 1})
        FreeModuleElement({'a': 1, 'b': 1})

        """
        answer = type(self)(other)
        answer.__dict__ = self.__dict__
        answer.preferred_rep()
        return answer

    def zero(self):
        """Initializes 0 with same type and attribute values as *self*.

        RETURNS
        -------
        type(*self*) object
            The initialized empty object

        EXAMPLE
        -------
        >>> x = FreeModuleElement({'a': 1})
        >>> x + x.zero() == x
        True

        """
        return self.create()
