from math import factorial

from ..barratt_eccles import BarrattEcclesElement
from ..cubical import Cubical
from ..simplicial import Simplicial
from ..surjection import SurjectionElement
from ..symmetric import SymmetricRing


class Steenrod:
    """Produces Steenrod operation elements and chains."""

    @staticmethod
    def _check_model(model):
        models = {'surjection', 'barratt_eccles'}
        if model not in models:
            raise ValueError(
                f"model must be one of {sorted(models)}, not {model!r}")

    @staticmethod
    def may_structure(arity, degree, model='surjection', torsion=None,
                      convention=None):
        r"""Representative of the requested May-Steenrod structure.

        PARAMETERS
        ----------
        arity : :class:`int`
            The arity considered.
        degree : :class:`int`
            The degree considered.
        model : :class:`str`, default ``'surjection'``
            Either ``'surjection'`` or ``'barratt_eccles'``.
        torsion : :class:`int` or ``None``, default ``None``
            The integer :math:`n` defining
            :math:`\mathbb Z/n \mathbb Z`, with :math:`n = 0` representing
            :math:`\mathbb Z`. Negative values are normalized to their
            absolute value.
        convention : :class:`str` or ``None``, default ``None``
            The sign convention of the surjection operad. Only used for
            ``model='surjection'``.

        RETURNS
        -------
        :class:`comch.surjection.SurjectionElement` or
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The image of the preferred basis element under the chosen
            May-Steenrod structure.

        """
        Steenrod._check_model(model)
        if model == 'surjection':
            return Steenrod._surjection_may_structure(
                arity, degree, torsion=torsion, convention=convention)
        if convention is not None:
            raise ValueError(
                "convention is only defined for model='surjection'")
        return Steenrod._barratt_eccles_may_structure(
            arity, degree, torsion=torsion)

    @staticmethod
    def operation(p, s, q, bockstein=False, model='surjection',
                  convention=None):
        r"""Chain level representative of a Steenrod operation.

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
        model : :class:`str`, default ``'surjection'``
            Either ``'surjection'`` or ``'barratt_eccles'``.
        convention : :class:`str` or ``None``, default ``None``
            The sign convention of the surjection operad. Only used for
            ``model='surjection'``. If ``None``, the McClure-Smith convention
            is used for the surjection model.

        RETURNS
        -------
        :class:`comch.surjection.SurjectionElement` or
        :class:`comch.barratt_eccles.BarrattEcclesElement`
            The chosen operad representative of the Steenrod operation.

        """
        Steenrod._check_model(model)
        if not all(isinstance(i, int) for i in {p, s, q}):
            raise TypeError('initialize with three int p,s,n')
        if not isinstance(bockstein, bool):
            raise TypeError('bockstein must be a boolean')
        if p == 2 and bockstein:
            raise TypeError('bP only defined for odd primes')
        if model == 'barratt_eccles' and convention is not None:
            raise ValueError(
                "convention is only defined for model='surjection'")

        if p == 2:
            coeff = 1
            degree = s - q
            if degree < 0:
                if model == 'surjection':
                    return SurjectionElement(torsion=2)
                return BarrattEcclesElement(torsion=2)
        else:
            b = int(bockstein)
            coeff = (-1) ** (q // 2 + s)
            if q % 2:
                coeff *= factorial((p - 1) // 2)
            degree = (2 * s - q) * (p - 1) - b
            if degree < 0:
                if model == 'surjection':
                    return SurjectionElement(torsion=p)
                return BarrattEcclesElement(torsion=p)

        if model == 'surjection':
            if convention is None:
                convention = 'McClure-Smith'
            element = Steenrod.may_structure(
                p, degree, model=model, torsion=p, convention=convention)
        else:
            element = Steenrod.may_structure(
                p, degree, model=model, torsion=p)
        return int(coeff) * element

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

        surj = Steenrod.operation(
            p, s, q, bockstein=bockstein, model='surjection',
            convention=convention)
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

    @staticmethod
    def _surjection_may_structure(arity, degree, torsion=None,
                                  convention=None):
        def i(surj, iterate=1):
            if iterate == 1:
                answer = surj.zero()
                for k, v in surj.items():
                    answer += answer.create(
                        {(1,) + tuple(j + 1 for j in k): v})
                return answer
            if iterate > 1:
                return i(i(surj, iterate=iterate - 1))

        def p(surj, iterate=1):
            if iterate == 1:
                answer = surj.zero()
                for k, v in surj.items():
                    if k.count(1) == 1:
                        idx = k.index(1)
                        new_k = (tuple(j - 1 for j in k[:idx]) +
                                 tuple(j - 1 for j in k[idx + 1:]))
                        answer += answer.create({new_k: v})
                return answer
            if iterate > 1:
                return p(p(surj, iterate=iterate - 1))

        def s(surj):
            answer = surj.zero()
            for k, v in surj.items():
                answer += answer.create({(1,) + tuple(j for j in k): v})
            return answer

        def h(surj):
            answer = s(surj)
            for r in range(1, arity - 1):
                answer += i(s(p(surj, r)), r)
            return answer

        operators = {
            0: SymmetricRing.norm_element(arity),
            1: SymmetricRing.transposition_element(arity)
        }

        def psi(arity, degree, convention=convention):
            if degree == 0:
                return SurjectionElement({tuple(range(1, arity + 1)): 1},
                                         convention=convention)
            previous = psi(arity, degree - 1, convention=convention)
            acted_on = operators[degree % 2] * previous
            return h(acted_on)

        if convention is None:
            convention = SurjectionElement.default_convention
        integral_answer = psi(arity, degree, convention=convention)
        if torsion is not None:
            integral_answer.set_torsion(torsion)
        return integral_answer

    @staticmethod
    def _barratt_eccles_may_structure(arity, degree, torsion=None):
        operators = {
            0: SymmetricRing.norm_element(arity),
            1: SymmetricRing.transposition_element(arity)
        }

        def psi(arity, degree):
            if degree == 0:
                return BarrattEcclesElement(
                    {(tuple(range(1, arity + 1)),): 1})
            previous = psi(arity, degree - 1)
            acted_on = operators[degree % 2] * previous
            identity = tuple(range(1, arity + 1))
            new_data = {(identity,) + k: v for k, v in acted_on.items()}
            return BarrattEcclesElement(new_data)

        integral_answer = psi(arity, degree)
        if torsion is not None:
            integral_answer.set_torsion(torsion)
        return integral_answer
