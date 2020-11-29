from clesto import *
from _utils import pairwise
from itertools import chain, combinations, product
from operator import itemgetter
from functools import reduce


class Surjection_element(Surjection_element):

    def __call__(self, other):
        '''Action on an basis element in the normalized chains of a standard
        cube or simplex represented by an arity 1 element in the (cubical)
        Eilenberg-Zilber operad.

        Examples
        --------

        # chain map check

        >>> s = Surjection_element({(3, 2, 1, 3, 1, 2): 1},\
                                       convention='McClure-Smith')
        >>> x = EilenbergZilber_element({((0, 1, 2, 3), ): 1})
        >>> ds_x = s.boundary()(x)
        >>> d_sx = s(x).boundary()
        >>> sdx = s(x.boundary())
        >>> d_sx - ((-1)**(s.degree)) * sdx == ds_x
        True

        '''

        def _sign(k1, k2):
            '''...
            '''
            def ordering_sign(permu, weights):
                '''Returns the exponent of the Koszul sign of the given
                permutation acting on the elements of degrees given by the
                list of weights

                '''
                sign_exp = 0
                for idx, j in enumerate(permu):
                    to_add = [weights[permu.index(i)] for
                              i in permu[idx + 1:] if i < j]
                    sign_exp += weights[idx] * sum(to_add)
                return sign_exp % 2

            def action_sign(ordered_k1, ordered_weights):
                '''Given a ordered tuple [1,..,1, 2,...,2, ..., r,...,r]
                and weights [w_1, w_2, ..., w_{r+d}] of the same length, gives
                the kozul sign obtained by inserting from the left a weight 1
                operator between equal consecutive elements.

                '''
                sign_exp = 0
                for idx, (i, j) in enumerate(pairwise(ordered_k1)):
                    if i == j:
                        sign_exp += sum(ordered_weights[:idx + 1])
                return sign_exp % 2

            sign_exp = 0
            weights = [e.dimension % 2 for e in k2]
            inv_ordering_permu = [pair[0] for pair in
                                  sorted(enumerate(k1), key=itemgetter(1))]
            ordering_permu = tuple(inv_ordering_permu.index(i)
                                   for i in range(len(inv_ordering_permu)))
            sign_exp += ordering_sign(ordering_permu, weights)
            ordered_k1 = list(sorted(k1))
            ordered_weights = [weights[i] for i in inv_ordering_permu]
            sign_exp += action_sign(ordered_k1, ordered_weights)
            return (-1)**sign_exp

        def _simplicial(self, other):
            answer = other.zero()
            pre_join = other.iterated_diagonal(self.arity + self.degree - 1)
            for (k1, v1), (k2, v2) in product(self.items(), pre_join.items()):
                new_k = []
                zero_summand = False
                for i in range(1, max(k1) + 1):
                    to_join = (spx for idx, spx in enumerate(k2)
                               if k1[idx] == i)
                    joined = Simplex(reduce(lambda x, y: x + y, to_join))
                    if joined.is_degenerate():
                        zero_summand = True
                        break
                    new_k.append(joined)

                if not zero_summand:
                    if self.torsion == 2:
                        sign = 1
                    else:
                        sign = _sign(k1, k2)
                    answer += answer.create({tuple(new_k): sign * v1 * v2})
            return answer

        def _cubical(self, other):
            pass

        if not self or not other:
            return other.zero()

        if other.arity != 1:
            raise TypeError(f'action only on arity 1, not {other.arity}')

        if self.degree is None or self.arity is None:
            raise TypeError('defined for homogeneous surjections')

        if self.torsion != other.torsion:
            raise TorsionError

        if isinstance(other, EilenbergZilber_element):
            if self.convention != 'McClure-Smith':
                raise NotImplementedError
            return _simplicial(self, other)

        elif isinstance(other, CubicalEilenbergZilber_element):
            return _cubical(self, other)

        else:
            raise NotImplementedError


s = Surjection_element({(1, 2, 1): 1}, convention='McClure-Smith')
ez = EilenbergZilber.standard_element(2)
print('answer', s(ez))
