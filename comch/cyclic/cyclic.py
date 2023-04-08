from ..free_module import FreeModuleElement


class CyclicGroupElement:
    r"""Element in a finite cyclic group.

    We refer to elements in the cyclic group of order :math:`r`, denoted
    :math:`\mathrm C_r`, as cyclic group elements.
    We identify :math:`\mathrm C_r` with :math:`\mathbb Z/r\mathbb Z`.

    """

    default_order = 1
    """Class attribute: Used if :attr:`order` is ``None`` during
    initialization."""

    def __init__(self, value=None, order=None):
        """Initializes *self*.

        PARAMETERS
        ----------
        value : :class:'int'
            An integer representing a class in :math:`mathrm C_r`.

        order : :class:'int'
            Positive integer representing the order :math:`r`.

        EXAMPLE
        -------
        >>> print(CyclicGroupElement(5, 3))
        2

        """
        if order is None:
            order = CyclicGroupElement.default_order
        self.value = value % order
        self.order = order

    def __str__(self):
        return str(self.value)

    def __mul__(self, other):
        r"""Product: *self* * *other*.

        This product agrees with the modular sum of integer
        *self* :math:`+` *other*  (mod r).

        RETURNS
        -------
        :class: `comch.cyclic.CyclicGroupElement` object
            The product of *self* and *other*.

        EXAMPLE
        -------
        >>> print(CyclicGroupElement(5, 3) * CyclicGroupElement(2, 3))
        1

        """
        if not isinstance(other, CyclicGroupElement) or self.order != other.order:
            raise ValueError("Both operands must be CyclicGroupElement instances "
                             "with the same order.")
        return CyclicGroupElement((self.value + other.value) % self.order, self.order)

    def __int__(self):
        return int(self.value)
class CyclicModuleElement(FreeModuleElement):
    r"""Element in :math:`N E\mathrm C_r`

    For a non-negative integer :math:`r` define the simplicial set
    :math:`E(\mathrm C_r)` by

    .. math::

       \begin{aligned}
       E(\mathrm C_r)_n &=
       \{ (\sigma_0, \dots, \sigma_n)\ |\ \sigma_i \in \mathrm C_r\}, \\
       d_i(\sigma_0, \dots, \sigma_n) &=
       (\sigma_0, \dots, \widehat{\sigma}_i, \dots, \sigma_n), \\
       s_i(\sigma_0, \dots, \sigma_n) &=
       (\sigma_0, \dots, \sigma_i, \sigma_i, \dots, \sigma_n),
       \end{aligned}

    corresponding to the unreduced bar construction on the monoid
    :math:`\mathrm C_r`. It is equipped with a left
    :math:`\mathrm C_r`-action defined on basis elements by

    .. math::
        \sigma (\sigma_0, \dots, \sigma_n) =
        (\sigma \sigma_0, \dots, \sigma \sigma_n).

    The chain complex resulting from applying the functor of integral
    normalized chains to it is denoted :math:`\mathcal C(r)`.

    """

    def __init__(self, data=None, order=None, torsion=None):
        """Initializes *self*.

        PARAMETERS
        ----------
        data : :class:`dict` or ``None``, default: ``None``
            Dictionary representing a linear combination of basis elements. Items
            in the dictionary correspond with pairs `basis_element: coefficient`.
            Each basis_element must create a :class:`tuple` of
            :class:`commch.cyclic.CyclicGroupElement` of the same order
            and `coefficient` must be an :class:`int`.
        order : :class:`int`
            The order :math:`r` of the group :math:`\mathrm C_r`.
        torsion : :class:`int`
            The non-neg int :math:`n` defining :math:`\mathbb Z/n \mathbb Z`.

        EXAMPLES
        --------
        >>> x = CyclicModuleElement()
        >>> print(x)
        0
        >>> y = CyclicModuleElement({(1, 2): 1, (2, 3): -1}, order=2, torsion=3)
        >>> print(y)
        (1,0) + 2(0,1)

        """
        if data is None and order is None:
            order = CyclicGroupElement.default_order

        if data:
            first_key = next(iter(data))
            are_group_elements = isinstance(first_key[0], CyclicGroupElement)

            # Determine the order if not given
            if order is None:
                if are_group_elements:
                    order = first_key[0].order
                else:
                    order = CyclicGroupElement.default_order

            # Make the input into dict of tuple of CGE
            if not are_group_elements:
                new_data = {}
                for k, v in data.items():
                    new_k = tuple(CyclicGroupElement(int(i), order) for i in k)
                    new_data[new_k] = v
                data = new_data

        self.order = order
        super(CyclicModuleElement, self).__init__(data=data, torsion=torsion)

    def __str__(self):
        new_data = {}
        for k, v in self.items():
            new_k = tuple(int(i) for i in k)
            new_data[new_k] = v
        string = str(FreeModuleElement(new_data, torsion=self.torsion))
        return string.replace(', ', ',')
