from collections import Counter


class Module_element(Counter):
    """Elements in a free module over :math:`\\mathbb{Z}` or
    :math:`\\mathbb{Z}/n \\mathbb{Z}`.

    Parameters
    ----------

    data : dict or None, default: ``None``
        Dictionary representing a linear cobination of basis elements.
        Items in the dict correspond with pairs (basis element: coefficient).
    torsion : int or 'free', default 'free'
        The torsion of the underlying ring :math:`\\mathbb{Z}` or
        :math:`\\mathbb{Z}/n \\mathbb{Z}`

    Attributes
    ----------

    torsion : int or 'free', default 'free'
        The torsion of the underlying ring :math:`\\mathbb{Z}` or
        :math:`\\mathbb{Z}/n \\mathbb{Z}`

    Example
    -------

    >>> print(Module_element())
    0
    >>> print(Module_element({'a': 1, 'b': -1, 'c': 0}))
    a - b

    """

    default_torsion = 'free'

    def __init__(self, data=None, torsion=None):

        if torsion is None:
            torsion = type(self).default_torsion

        self.torsion = torsion

        super(Module_element, self).__init__(data)

        self._reduce_rep()

    def __hash__(self):
        return hash(frozenset(self))

    def __str__(self):
        """Coefficient first representation."""
        if not self:
            return '0'
        else:
            answer = ''
            for key, value in self.items():
                if value < -1:
                    answer += f'- {abs(value)}{key} '
                elif value == -1:
                    answer += f'- {key} '
                elif value == 1:
                    answer += f'+ {key} '
                elif value > 1:
                    answer += f'+ {value}{key} '
            if answer[0] == '+':
                answer = answer[2:]

            return answer[:-1]

    def __add__(self, other):
        """Addition: *self* + *other*.

        Parameters
        ----------

        other : :class:`clesto.basics.module.Module_element` object
            The element to add to *self*.

        Returns
        -------

        :class:`clesto.basics.module.Module_element` object
            The sum of *self* and *other*.

        Example
        -------

        >>> Module_element({'a': 1, 'b': 2}) + Module_element({'a': 1})
        Module_element({'a': 2, 'b': 2})

        """
        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')
        answer = self.create(self)
        answer.update(other)
        answer._reduce_rep()
        return answer

    def __sub__(self, other):
        """Diference: *self* - *other*.

        Parameters
        ----------

        other : :class:`clesto.basics.module.Module_element` object
            The element to substract from *self*.

        Returns
        -------

        :class:`clesto.basics.module.Module_element` object
            The difference of *self* and *other*.

        Example
        -------

        >>> Module_element({'a': 1, 'b': 2}) - Module_element({'a': 1})
        Module_element({'b': 2})

        """
        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')
        answer = self.create(self)
        answer.subtract(other)
        answer._reduce_rep()
        return answer

    def __rmul__(self, c):
        """Scaling: *c* * *self*.

        Parameters
        ----------

        other : int
            The element to scale *self*. by.

        Returns
        -------

        :class:`clesto.basics.module.Module_element` object
            The scaling of *self* by *other*.

        Example
        -------

        >>> 3 * Module_element({'a':1, 'b':2})
        Module_element({'b': 6, 'a': 3})

        """
        if not isinstance(c, int):
            raise TypeError(f'Act only by int not by type {type(c)}')

        scaled = {k: c * v for k, v in self.items()}
        return self.create(scaled)

    def __neg__(self):
        """Additive inverse: - *self*.

        Returns
        -------

        :class:`clesto.basics.module.Module_element` object
            the additive inverse of *self*.

        Example
        -------

        >>> - Module_element({'a': 1, 'b': 2})
        Module_element({'a': -1, 'b': -2})

        """
        return self.__rmul__(-1)

    def __iadd__(self, other):
        """In place addition: *self* += *other*.

        Parameters
        ----------

        other : :class:`clesto.basics.module.Module_element` object
            The element to add to *self*.

        Example
        -------

        >>> x = Module_element({'a': 1, 'b': 2})
        >>> x += Module_element({'a': 3, 'b': 6})
        >>> x
        Module_element({'b': 8, 'a': 4})

        """
        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')
        self.update(other)
        self._reduce_rep()
        return self

    def __isub__(self, other):
        """In place difference: *self* -= *other*.

        Parameters
        ----------

        other : :class:`clesto.basics.module.Module_element` object
            The element to substract from *self*.

        Example
        -------

        >>> x = Module_element({'a': 1, 'b': 2})
        >>> x -= Module_element({'a': 3, 'b': 6})
        >>> x
        Module_element({'a': -2, 'b': -4})

        """
        if self.torsion != other.torsion:
            raise TypeError('only defined for equal attribute torsion')
        self.subtract(other)
        self._reduce_rep()
        return self

    def _reduce_rep(self):
        """The preferred representative of *self*.

        The preferred representative has coefficient in {0,...,torsion-1}
        if attribute torsion is not 'free', and no pairs key:value with
        value = 0.

        Example
        -------

        >>> Module_element({'a': 1, 'b': 2, 'c': 0})
        Module_element({'b': 2, 'a': 1})

        """
        # reducing coefficients mod torsion
        if self.torsion != 'free':
            for key, value in self.items():
                self[key] = value % self.torsion

        # removing key:value pairs with value = 0
        zeros = [k for k, v in self.items() if not v]
        for key in zeros:
            del self[key]

    def set_torsion(self, torsion):
        """Sets the torsion of *self*.

        Parameters
        ----------

        torsion : int or 'free'
            The new `torsion` of *self*

        Example
        -------

        >>> Module_element({'a': 1, 'b': 2}).set_torsion(2)
        Module_element({'a': 1})

        """
        setattr(self, 'torsion', torsion)
        self._reduce_rep()
        return self

    def create(self, other=None):
        """Instantiates data with same type and attribute values as *self*.

        Parameters
        ----------

        other : dict or None, default: ``None``
            Data to be initialized.

        Returns
        -------

        created : type(*self*) object
            The initialized object with the given data

        Example
        -------

        >>> x =  Module_element({'a': 1})
        >>> x + x.create({'b': 1})
        Module_element({'a': 1, 'b': 1})

        """
        answer = type(self)(other)
        answer.__dict__ = self.__dict__
        return answer

    def zero(self):
        """Instantiates 0 with same type and attribute values as *self*.

        Returns
        -------

        created : type(*self*) object
            The initialized empty object

        Example
        -------

        >>> x = Module_element({'a': 1})
        >>> x + x.zero() == x
        True

        """
        return self.create()
