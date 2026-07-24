from ..barratt_eccles import BarrattEccles
from ..cubical import Cubical
from ..simplicial import Simplicial
from ..surjection import Surjection


class Steenrod:
    r"""Chain-level representatives related to Steenrod operations.

    This namespace collects convenience methods for evaluating universal
    chain representatives of Steenrod operations on standard simplices and
    cubes, and Cartan and Adem coboundaries on standard simplices.

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
    def cartan_coboundaries(
            i, n, shape='simplex', convention=None, arity=2):
        r"""Chain representatives of Cartan coboundaries.

        At the even prime, the Cartan formula for Steenrod squares is
        equivalent to the existence of natural cochains
        :math:`\zeta_i(\alpha \otimes \beta)` satisfying

        .. math::

           \delta \zeta_i(\alpha \otimes \beta) =
           (\alpha \smallsmile_0 \beta) \smallsmile_i
           (\alpha \smallsmile_0 \beta)
           + \sum_{i=j+k}
           (\alpha \smallsmile_j \alpha) \smallsmile_0
           (\beta \smallsmile_k \beta).

        Such a cochain is called a Cartan :math:`i`-coboundary. At an odd
        prime :math:`p`, the effective relators of [CMM] give the cochain
        homotopies underlying

        .. math::

           P_s(xy) = \sum_{a+b=s}P_a(x)P_b(y)

        and

        .. math::

           \beta P_s(xy) =
           \sum_{a+b=s}
           \left(\beta P_a(x)P_b(y)
           + (-1)^{|x|}P_a(x)\beta P_b(y)\right).

        This method returns the value on the top :math:`n`-simplex of the
        chain map :math:`z_i \colon C \to C^{\otimes 2r}` obtained by acting
        with the arity-:math:`r` Cartan homotopy. For ``arity=2`` it is
        characterized by

        .. math::

           \zeta_i(\alpha \otimes \beta)(c) =
           (\alpha \otimes \alpha \otimes \beta \otimes \beta)z_i(c).

        On simplicial chains the action factors through the table reduction
        morphism to the surjection operad. The Barratt-Eccles action on
        cubical chains is not yet implemented.

        Since :math:`\alpha` and :math:`\beta` can be assumed homogeneous, the
        output is filtered to the summands where the first :math:`r` tensor
        factors have equal degree and the last :math:`r` tensor factors have
        equal degree.

        PARAMETERS
        ----------
        i : :class:`int`
            The index of the Cartan coboundary.
        n : :class:`int`
            The degree of the universal basis element.
        shape : :class:`str`, default ``'simplex'``
            Either ``'simplex'`` or ``'cube'``. The latter currently raises
            :class:`NotImplementedError`.
        convention : :class:`str` or ``None``, default ``None``
            The sign convention of the surjection operad, used when *shape*
            is ``'simplex'``. If ``None``, the McClure-Smith convention is
            used.
        arity : :class:`int`, default ``2``
            The arity of the cyclic resolution. Use an odd prime to obtain
            the corresponding odd-prime Cartan coboundary.

        RETURNS
        -------
        :class:`comch.simplicial.SimplicialElement`
            The tensor product element representing the Cartan
            :math:`i`-coboundary on the selected universal basis element.

        REFERENCES
        ----------
        [Med20]: A. M. Medina-Mardones. "An effective proof of the Cartan
        formula: the even prime." Journal of Pure and Applied Algebra 224
        (2020), no. 12, 106444.
        https://doi.org/10.1016/j.jpaa.2020.106444

        [CMM]: F. Cantero-Morán and A. M. Medina-Mardones. "An effective
        proof of the Cartan formula: odd primes." Homology, Homotopy and
        Applications 27 (2025), no. 1, 207--234.
        https://doi.org/10.4310/HHA.2025.v27.n1.a12

        """

        if convention is None:
            convention = (
                'McClure-Smith' if arity == 2 else 'Berger-Fresse')

        if shape == 'simplex':
            surj = Surjection.cartan_homotopy(
                i, convention=convention, arity=arity)
            element = Simplicial.standard_element(n, torsion=arity)
            answer = surj(element)
        elif shape == 'cube':
            homotopy = BarrattEccles.cartan_homotopy(i, arity=arity)
            element = Cubical.standard_element(n, torsion=arity)
            answer = homotopy(element)
        else:
            raise ValueError("shape must be 'simplex' or 'cube'")

        return Steenrod._filter_bihomogeneous(answer, arity=arity)

    @staticmethod
    def adem_coboundaries(q, p, n, shape='simplex', convention=None):
        r"""Chain representatives of Adem coboundaries at the even prime.

        For integers :math:`q > p \geq 0`, this method returns the value on
        the top :math:`n`-simplex of the chain map
        :math:`a_{q,p} \colon C \to C^{\otimes 4}` obtained by acting with
        the complete Adem coboundary. It includes both the Adem homotopy and
        the two non-symmetric corrections of [BMM, Theorem 0.2.5], and is
        characterized by

        .. math::

           A_{q,p}(\alpha)(c) =
           (\alpha \otimes \alpha \otimes \alpha \otimes \alpha)
           a_{q,p}(c).

        The parameters :math:`q,p` index the universal classes used in the
        chain-level construction; they are not, in general, the exponents
        of a single composite :math:`\mathrm{Sq}^q\mathrm{Sq}^p`.

        On simplicial chains the action factors through table reduction to
        the surjection operad. The Barratt-Eccles action on cubical chains
        is not yet implemented.

        Since :math:`\alpha` can be assumed homogeneous, the output is
        filtered to the summands where all four tensor factors have equal
        degree.

        PARAMETERS
        ----------
        q, p : :class:`int`
            Non-negative integers with ``q > p``.
        n : :class:`int`
            The degree of the universal basis element.
        shape : :class:`str`, default ``'simplex'``
            Either ``'simplex'`` or ``'cube'``. The latter currently raises
            :class:`NotImplementedError`.
        convention : :class:`str` or ``None``, default ``None``
            The sign convention of the surjection operad, used when *shape*
            is ``'simplex'``. If ``None``, the McClure-Smith convention is
            used.

        RETURNS
        -------
        :class:`comch.simplicial.SimplicialElement` or
        :class:`comch.cubical.CubicalElement`
            The tensor product element representing the Adem
            :math:`(q,p)`-coboundary on the selected universal basis element.

        """

        if convention is None:
            convention = 'McClure-Smith'

        if shape == 'simplex':
            surj = Surjection.adem_coboundary(q, p, convention=convention)
            element = Simplicial.standard_element(n, torsion=2)
            answer = surj(element)
        elif shape == 'cube':
            coboundary = BarrattEccles.adem_coboundary(q, p)
            element = Cubical.standard_element(n, torsion=2)
            answer = coboundary(element)
        else:
            raise ValueError("shape must be 'simplex' or 'cube'")

        return Steenrod._filter_homogeneous(answer)

    @staticmethod
    def _filter_homogeneous(element):
        homogeneous = {}
        for k, v in element.items():
            if len(set(elmt.dimension for elmt in k)) == 1:
                homogeneous[k] = v
        return element.create(homogeneous)

    @staticmethod
    def _filter_bihomogeneous(element, arity=2):
        homogeneous = {}
        for k, v in element.items():
            if len(k) != 2 * arity:
                continue
            left_dimensions = {elmt.dimension for elmt in k[:arity]}
            right_dimensions = {elmt.dimension for elmt in k[arity:]}
            if len(left_dimensions) == len(right_dimensions) == 1:
                homogeneous[k] = v
        return element.create(homogeneous)
