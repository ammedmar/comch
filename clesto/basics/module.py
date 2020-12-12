from collections import Counter


class Module_element(Counter):
    """Elements in a free module over the (modular) integers.

    """

    default_torsion = 'free'

    def __init__(self, data=None, torsion=None):
        """Initialize an instance of Module_element

        Create a new, empty Module_element object representing 0, and, if given,
        initialize a Module_element from a dict with integer values.

        >>> print(Module_element())
        0
        >>> print(Module_element({'a': 1, 'b': -1, 'c': 0}))
        a - b

        """
        if torsion is None:
            torsion = type(self).default_torsion
        self.torsion = torsion

        super(Module_element, self).__init__(data)

        self._reduce_rep()

    def __str__(self):
        """Coefficient first representation.

        >>> str(Module_element({'a': 1, 'b': -1}))
        'a - b'

        """
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

    def __hash__(self):
        return hash(frozenset(self))

    def __add__(self, other):
        """Addition: self + other

        >>> Module_element({'a': 1, 'b': 2}) + Module_element({'a': 1})
        Module_element({'a': 2, 'b': 2})

        """
        if self.torsion != other.torsion:
            raise TypeError('Unequal torsion attribute')
        answer = self.create(self)
        answer.update(other)
        answer._reduce_rep()
        return answer

    def __sub__(self, other):
        """Diference: self - other

        >>> Module_element({'a': 1, 'b': 2}) - Module_element({'a': 1})
        Module_element({'b': 2})

        """
        if self.torsion != other.torsion:
            raise TypeError('Unequal torsion attribute')
        answer = self.create(self)
        answer.subtract(other)
        answer._reduce_rep()
        return answer

    def __rmul__(self, c):
        """Integral scaling: c * self

        >>> 3 * Module_element({'a':1, 'b':2})
        Module_element({'b': 6, 'a': 3})

        """
        if not isinstance(c, int):
            raise TypeError(f'can not act by non-int of type {type(c)}')

        scaled = {k: c * v for k, v in self.items()}
        return self.create(scaled)

    def __neg__(self):
        """Additive inverse: - self

        >>> - Module_element({'a': 1, 'b': 2})
        Module_element({'a': -1, 'b': -2})

        """
        return self.__rmul__(-1)

    def __iadd__(self, other):
        """In place addition: self + other

        >>> x = Module_element({'a': 1, 'b': 2})
        >>> x += Module_element({'a': 3, 'b': 6})
        >>> x
        Module_element({'b': 8, 'a': 4})

        """
        if self.torsion != other.torsion:
            raise TypeError('Unequal torsion attribute')
        self.update(other)
        self._reduce_rep()
        return self

    def __isub__(self, other):
        """In place difference: self - other

        >>> x = Module_element({'a': 1, 'b': 2})
        >>> x -= Module_element({'a': 3, 'b': 6})
        >>> x
        Module_element({'a': -2, 'b': -4})

        """
        if self.torsion != other.torsion:
            raise TypeError('Unequal torsion attribute')
        self.subtract(other)
        self._reduce_rep()
        return self

    def _reduce_rep(self):
        """The preferred representative of the free module element

        It reduces all values mod n if torsion is n and removes
        key:value pairs with value = 0.

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
        """Sets the torsion of an element.

        >>> Module_element({'a': 1, 'b': 2}).set_torsion(2)
        Module_element({'a': 1})

        """
        setattr(self, 'torsion', torsion)
        self._reduce_rep()
        return self

    def create(self, other=None):
        """Instantiates data with the same type and attribute values as self.

        >>> x =  Module_element({'a': 1})
        >>> x + x.create({'b': 1})
        Module_element({'a': 1, 'b': 1})

        """
        answer = type(self)(other)
        answer.__dict__ = self.__dict__
        return answer

    def zero(self):
        """Instantiates 0 with same type and attribute values as self.

        >>> x = Module_element({'a': 1})
        >>> x + x.zero() == x
        True

        """
        return self.create()
