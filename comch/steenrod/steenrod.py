from ..cubical import Cubical
from ..simplicial import Simplicial
from ..surjection import Surjection


class Steenrod:
    r"""Produces Steenrod operation.

    For convenience, we recall the construction. Let :math:`A` be an algebra
    over an operad equipped with a May-Steenrod map, so there is a composite

    .. math::
       \psi : W \longrightarrow \mathcal O \longrightarrow \mathrm{End}_A.

    For a prime :math:`p`, define the linear map
    :math:`D_d : A \otimes \mathbb F_p \to A \otimes \mathbb F_p` by

    .. math::
       D_d(a) = \begin{cases}
       \psi(e_d)(a^{\otimes p})& d \geq 0 \\
       0 & d < 0.
       \end{cases}

    For any integer :math:`s`, the Steenrod operations

    .. math::
        P_s : H_\bullet(A; \mathbb F_2)
        \to
        H_{\bullet + s}(A; \mathbb F_2)

    and, for :math:`p > 2`,

    .. math::
       P_s : H_\bullet(A; \mathbb F_p)
       \to
       H_{\bullet + 2s(p-1)}(A; \mathbb F_p)
       \qquad
       \beta P_s : H_\bullet(A; \mathbb F_p)
       \to
       H_{\bullet+2s(p-1)-1}(A; \mathbb F_p)

    are defined for a class :math:`[a]` of degree
    :math:`q` respectively by

    .. math:: P_s\big([a]\big) = \big[D_{s-q}(a)\big] \qquad

    and

    .. math::
       P_s\big([a]\big) = \big[(-1)^s \nu(q) D_{(2s-q)(p-1)}(a)\big]
       \qquad
       \beta P_s\big([a]\big) = \big[(-1)^s \nu(q)D_{(2s-q)(p-1)-1}(a)\big]

    where :math:`\nu(q) = (-1)^{q(q-1)m/2}(m!)^q` and :math:`m = (p-1)/2`.

    """

    @staticmethod
    def chain(p, s, q, bockstein=False, shape='simplex', convention=None):
        """Chain representative of a Steenrod operation.

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

        def filter_homogeneous(element):
            homogeneous = {}
            for k, v in element.items():
                if len(set(elmt.dimension for elmt in k)) == 1:
                    homogeneous[k] = v
            return element.create(homogeneous)

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

        return filter_homogeneous(surj(element))
