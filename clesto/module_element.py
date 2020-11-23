from collections import Counter


class TorsionError(Exception):
    """Exception raised for unequal torsion.

    Attributes:
        message : explanation of the error
    """

    def __init__(self, message='attribute torsion must be the same'):
        self.message = message
        super(TorsionError, self).__init__(message)


class Module_element(Counter):
    '''
    Counter with arithmetic improvements to handle (modular) integer values.

    Class constructed to model free module elements over the ring Z or Z/nZ.

    Parameters
    ----------
    data : dict

    Attributes
    ----------
    torsion : non-negative int or string 'free'.

    '''

    default_torsion = 'free'

    def __init__(self, data=None, torsion=None):
        '''...'''

        if torsion is None:
            torsion = type(self).default_torsion
        self.torsion = torsion

        super(Module_element, self).__init__(data)

        self._reduce_rep()

    def __str__(self):
        '''...'''

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
        '''The sum of two free module elements.

        >>> Module_element({'a': 1, 'b': 2}) + Module_element({'a': 1})
        Module_element({'a': 2, 'b': 2})

        '''
        if self.torsion != other.torsion:
            raise TorsionError
        self.update(other)
        return type(self)(self)

    def __sub__(self, other):
        '''The substraction of two free module elements.

        >>> Module_element({'a': 1, 'b': 2}) - Module_element({'a': 1})
        Module_element({'b': 2})

        '''
        if self.torsion != other.torsion:
            raise TorsionError
        self.subtract(other)
        return type(self)(self)

    def __rmul__(self, c):
        '''The scaling by c of a free module element.

        >>> 3*Module_element({'a':1, 'b':2})
        Module_element({'b': 6, 'a': 3})

        '''
        if not isinstance(c, int):
            raise TypeError(f'can not act by non-int of type {type(c)}')

        scaled = {k: c * v for k, v in self.items()}
        return self.create(scaled)

    def __neg__(self):
        '''The additive inverse of a free module element.

        >>> -Module_element({'a': 1, 'b': 2})
        Module_element({'a': -1, 'b': -2})

        '''
        return self.__rmul__(-1)

    def __iadd__(self, other):
        '''The in place addition of two free module elements.

        >>> x = Module_element({'a': 1, 'b': 2})
        >>> x += Module_element({'a': 3, 'b': 6})
        >>> x
        Module_element({'b': 8, 'a': 4})

        '''
        if self.torsion != other.torsion:
            raise TorsionError
        self.update(other)
        self._reduce_rep()
        return self

    def __isub__(self, other):
        '''The in place addition of two free module elements.

        >>> x = Module_element({'a': 1, 'b': 2})
        >>> x -= Module_element({'a': 3, 'b': 6})
        >>> x
        Module_element({'a': -2, 'b': -4})

        '''
        if self.torsion != other.torsion:
            raise TorsionError
        self.subtract(other)
        self._reduce_rep()
        return self

    def _reduce_rep(self):
        '''The preferred representative of the free module element.

        It reduces all values mod n if torsion is n and removes
        key:value pairs with value = 0.

        >>> Module_element({'a': 1, 'b': 2, 'c': 0})
        Module_element({'b': 2, 'a': 1})

        '''
        # reducing coefficients mod torsion
        if self.torsion != 'free':
            for key, value in self.items():
                self[key] = value % self.torsion

        # removing key:value pairs with value = 0
        zeros = [k for k, v in self.items() if not v]
        for key in zeros:
            del self[key]

    def set_torsion(self, torsion):
        '''...'''
        setattr(self, 'torsion', torsion)
        self._reduce_rep()
        return self

    def create(self, other=None):
        '''...'''
        return type(self)(other, torsion=self.torsion)

    def zero(self):
        '''...'''
        return self.create()

    def summands(self):
        '''...'''
        for k, v in self.items():
            yield self.create({k: v})


if __name__ == "__main__":
    import doctest
    doctest.testmod()
