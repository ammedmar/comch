from ..cubical import Cubical
from ..simplicial import Simplicial
from ..surjection import Surjection


class Steenrod:
    r"""Chain-level representatives related to Steenrod operations.

    This namespace collects convenience methods for evaluating universal
    chain representatives of Steenrod operations and their Cartan
    coboundaries on standard simplices and cubes.

    """

    @staticmethod
    def chain_operations(
            p, s, q, bockstein=False, shape='simplex', convention=None):
        r"""Chain representatives of Steenrod operations.

        Let :math:`A` be an algebra over an operad equipped with a
        May-Steenrod map, so that there is a composite

        .. math::

           \psi : W \longrightarrow \mathcal O \longrightarrow \mathrm{End}_A.

        For a prime :math:`p`, define

        .. math::

           D_d(a) = \begin{cases}
           \psi(e_d)(a^{\otimes p}) & d \geq 0, \\
           0 & d < 0.
           \end{cases}

        For any integer :math:`s`, the Steenrod operations

        .. math::

            P_s : H_q(A; \mathbb F_2)
            \to
            H_{q+s}(A; \mathbb F_2)

        and, for :math:`p > 2`,

        .. math::

           P_s : H_q(A; \mathbb F_p)
           \to
           H_{q + 2s(p-1)}(A; \mathbb F_p),
           \qquad
           \beta P_s : H_q(A; \mathbb F_p)
           \to
           H_{q+2s(p-1)-1}(A; \mathbb F_p)

        are defined, for a class :math:`[a]` of degree :math:`q`, by

        .. math::

           P_s\big([a]\big) = \big[D_{s-q}(a)\big]

        at :math:`p=2`, and by

        .. math::

           P_s\big([a]\big) =
           \big[(-1)^s \nu(q) D_{(2s-q)(p-1)}(a)\big],
           \qquad
           \beta P_s\big([a]\big) =
           \big[(-1)^s \nu(q)D_{(2s-q)(p-1)-1}(a)\big],

        where :math:`\nu(q) = (-1)^{q(q-1)m/2}(m!)^q` and
        :math:`m = (p-1)/2`.

        Simplicial and cubical cochains are examples of such algebras, with
        the structure induced by an :math:`E_\infty`-coalgebra on chains. By
        naturality, the cochain operations are controlled by their values on
        the universal basis: the top-degree generator of the chains on
        standard simplices or cubes.

        This method returns the value on that universal basis element of a
        chain map :math:`\mathbb P_s \colon C \to C^{\otimes p}` such that

        .. math::

           P_s[\alpha] = [\alpha^{\otimes p} \mathbb P_s(-)].

        Since :math:`\alpha` can be assumed homogeneous, the output is
        filtered to the homogeneous summands.

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
        shape : :class:`str`, default ``'simplex'``
            Either ``'simplex'`` or ``'cube'``.
        convention : :class:`str` or ``None``, default ``None``
            The sign convention of the surjection operad. If ``None``, the
            McClure-Smith convention is used.

        RETURNS
        -------
        :class:`comch.simplicial.SimplicialElement` or
        :class:`comch.cubical.CubicalElement`
            The tensor product element obtained by applying the representative
            of the specified Steenrod operation.

        """

        surj = Surjection.steenrod_operation_representative(
            p, s, q, bockstein=bockstein, convention=convention)
        b = int(bockstein)

        if p == 2:
            degree = q + s
        else:
            degree = q + 2 * s * (p - 1) - b

        if shape == 'simplex':
            element = Simplicial.standard_element(-degree, torsion=p)
        elif shape == 'cube':
            element = Cubical.standard_element(-degree, torsion=p)
        else:
            raise ValueError("shape must be 'simplex' or 'cube'")

        return Steenrod._filter_homogeneous(surj(element))

    @staticmethod
    def cartan_coboundaries(i, n, shape='simplex', convention=None):
        r"""Chain representatives of Cartan coboundaries at the even prime.

        The Cartan formula for Steenrod squares is equivalent to the existence
        of natural cochains :math:`\zeta_i(\alpha \otimes \beta)` satisfying

        .. math::

           \delta \zeta_i(\alpha \otimes \beta) =
           (\alpha \smallsmile_0 \beta) \smallsmile_i
           (\alpha \smallsmile_0 \beta)
           + \sum_{i=j+k}
           (\alpha \smallsmile_j \alpha) \smallsmile_0
           (\beta \smallsmile_k \beta).

        Such a cochain is called a Cartan :math:`i`-coboundary. This method
        returns the value on the top :math:`n`-simplex or :math:`n`-cube of the
        chain map :math:`z_i \colon C \to C^{\otimes 4}` obtained from the
        Cartan homotopy in the Barratt-Eccles operad followed by table
        reduction to the surjection operad. It is characterized by

        .. math::

           \zeta_i(\alpha \otimes \beta)(c) =
           (\alpha \otimes \alpha \otimes \beta \otimes \beta)z_i(c).

        Since :math:`\alpha` and :math:`\beta` can be assumed homogeneous, the
        output is filtered to the summands where the first two tensor factors
        have equal degree and the last two tensor factors have equal
        degree.

        PARAMETERS
        ----------
        i : :class:`int`
            The index of the Cartan coboundary.
        n : :class:`int`
            The degree of the universal basis element.
        shape : :class:`str`, default ``'simplex'``
            Either ``'simplex'`` or ``'cube'``.
        convention : :class:`str` or ``None``, default ``None``
            The sign convention of the surjection operad. If ``None``, the
            McClure-Smith convention is used.

        RETURNS
        -------
        :class:`comch.simplicial.SimplicialElement` or
        :class:`comch.cubical.CubicalElement`
            The tensor product element representing the Cartan
            :math:`i`-coboundary on the selected universal basis element.

        """

        if convention is None:
            convention = 'McClure-Smith'

        surj = Surjection.cartan_homotopy(i, convention=convention)

        if shape == 'simplex':
            element = Simplicial.standard_element(n, torsion=2)
        elif shape == 'cube':
            element = Cubical.standard_element(n, torsion=2)
        else:
            raise ValueError("shape must be 'simplex' or 'cube'")

        return Steenrod._filter_bihomogeneous(surj(element))

    @staticmethod
    def _filter_homogeneous(element):
        homogeneous = {}
        for k, v in element.items():
            if len(set(elmt.dimension for elmt in k)) == 1:
                homogeneous[k] = v
        return element.create(homogeneous)

    @staticmethod
    def _filter_bihomogeneous(element):
        homogeneous = {}
        for k, v in element.items():
            if k[0].dimension == k[1].dimension and \
                    k[2].dimension == k[3].dimension:
                homogeneous[k] = v
        return element.create(homogeneous)
