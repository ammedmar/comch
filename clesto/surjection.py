from clesto.module_element import Module_element, TorsionError
from clesto.symmetric import SymmetricGroup_element, \
    SymmetricModule_element, SymmetricModule, ArityError
from clesto.utils import pairwise, decompositions, distinct_permutations
from itertools import chain, combinations, product
from operator import attrgetter


class Surjection_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None, convention='Berger-Fresse'):
        '''...'''
        # check input data: dict with tuple of int keys
        if data:
            if not (isinstance(data, dict)
                    and all(isinstance(surj, tuple) for surj in data.keys())
                    and all(isinstance(i, int) for i in
                            chain.from_iterable(data.keys()))
                    ):
                raise TypeError(
                    'data type must be dict with tuple of int keys')

        if convention not in {'Berger-Fresse', 'McClure-Smith'}:
            raise TypeError(
                'convention must be either Berger-Fresse or McClure-Smith')

        # initialize element
        self.convention = convention

        super(Surjection_element, self).__init__(data=data, torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def zero(self):
        '''...'''
        return type(self)(torsion=self.torsion,
                          convention=self.convention)

    def create(self, other):
        '''...'''
        return type(self)(other, torsion=self.torsion,
                          convention=self.convention)

    @property
    def arity(self):
        '''...'''
        if not self:
            return None
        arities = set(max(surj) for surj in self.keys())
        if len(arities) == 1:
            return arities.pop()
        return None

    @property
    def degree(self):
        '''...'''
        if not self:
            return None
        degs = set(len(surj) - max(surj) for surj in self.keys())
        if len(degs) == 1:
            return degs.pop()
        return None

    @property
    def complexity(self):
        '''returns the complexity of an element in the Surjection operad

        >>> Surjection_element({(1, 2, 1, 3, 1, 3, 2, 3): 1}).complexity
        3
        >>> Surjection_element({(1, 2, 1):1, (1, 2): 1}).complexity
        2

        '''
        complexities = [0]
        for surjection in self.keys():
            for i, j in combinations(range(1, max(surjection) + 1), 2):
                r = tuple(k for k in surjection if k == i or k == j)
                cpxty = len([p for p, q in pairwise(r) if p != q])
                complexities.append(cpxty)

        return max(complexities)

    def boundary(self):
        '''boundary of self

        >>> s = Surjection_element({(1,2,1,3,1,3,2,3): 1}, torsion=2)
        >>> print(s.boundary())
        (2,1,3,1,3,2,3) + (1,2,3,1,3,2,3) + (1,2,1,3,1,2,3) + (1,2,1,3,1,3,2)
        >>> print(s.boundary().boundary())
        0

        >>> s = Surjection_element({(1,2,1,3,1,3,2,3): 1})
        >>> print(s.boundary())
        (2,1,3,1,3,2,3) + (1,2,3,1,3,2,3) + (1,2,1,3,1,2,3) - (1,2,1,3,1,3,2)
        >>> print(s.boundary().boundary())
        0

        >>> s = Surjection_element({(1,2,1,3,1,3,2,3): 1},
        ...                        convention='McClure-Smith')
        >>> print(s.boundary())
        (2,1,3,1,3,2,3) - (1,2,3,1,3,2,3) + (1,2,1,3,1,2,3) - (1,2,1,3,1,3,2)
        >>> print(s.boundary().boundary())
        0

        '''
        answer = self.zero()

        if self.torsion == 2:
            for k in self.keys():
                for idx in range(0, len(k)):
                    bdry_summand = k[:idx] + k[idx + 1:]
                    if k[idx] in bdry_summand:
                        answer += self.create({bdry_summand: 1})
            return answer

        if self.convention == 'Berger-Fresse':
            for k, v in self.items():
                # determining the signs of the summands
                signs = {}
                alternating_sign = 1
                for idx, i in enumerate(k):
                    if i in k[idx + 1:]:
                        signs[idx] = alternating_sign
                        alternating_sign *= (-1)
                    elif i in k[:idx]:
                        occurs = (pos for pos, j in enumerate(k[:idx]) if i == j)
                        signs[idx] = signs[max(occurs)] * (-1)
                    else:
                        signs[idx] = 0

                # computing the summands
                for idx in range(0, len(k)):
                    bdry_summand = k[:idx] + k[idx + 1:]
                    if k[idx] in bdry_summand:
                        answer += self.create({bdry_summand: signs[idx] * v})

        if self.convention == 'McClure-Smith':
            for k, v in self.items():
                sign = 1
                for i in range(1, max(k) + 1):
                    for idx in (idx for idx, j in enumerate(k) if j == i):
                        new_k = k[:idx] + k[idx + 1:]
                        if k[idx] in new_k:
                            answer += answer.create({new_k: v * sign})
                        sign *= -1
                    sign *= -1

        return answer

    def __rmul__(self, other):
        '''Left action by the appropriate symmetric group ring.

        # chain map checks:

        >>> rho = SymmetricModule.rotation_element(3)
        >>> surj = Surjection_element({(1, 2, 3, 1, 2): 1}, \
                                       convention='Berger-Fresse')
        >>> x, y = (rho * surj).boundary(), rho * surj.boundary()
        >>> x == y
        True

        >>> rho = SymmetricModule.rotation_element(3)
        >>> surj = Surjection_element({(1, 2, 3, 1, 3): 1}, \
                                       convention='McClure-Smith')
        >>> x, y = (rho * surj).boundary(), rho * surj.boundary()
        >>> x == y
        True

        '''

        def sign(perm, surj, convention):
            if convention == 'Berger-Fresse':
                return 1
            assert convention == 'McClure-Smith'
            signs = {0: 1, 1: -1}
            weights = [surj.count(i) - 1 for
                       i in range(1, max(surj) + 1)]
            answer = 0
            for idx, i in enumerate(perm):
                right = [weights[perm.index(j)] for
                         j in perm[idx + 1:] if i > j]
                answer += sum(right) * weights[idx]
            return signs[answer % 2]

        if isinstance(other, int):
            return super().__rmul__(other)

        if not isinstance(other, SymmetricModule_element):
            raise NotImplementedError

        if self.torsion != other.torsion:
            raise TorsionError

        if self.arity != other.arity:
            raise ArityError

        answer = self.zero()
        for k1, v1 in self.items():
            for k2, v2 in other.items():
                new_key = tuple(k2[i - 1] for i in k1)
                new_sign = sign(k2, k1, self.convention)
                answer += self.create({new_key: new_sign * v1 * v2})
        return answer

    def orbit(self, representation='trivial'):
        ''' Returns the preferred element in the symmetric orbit of an element

        >>> s = Surjection_element({(1, 3, 2): 1})
        >>> print(s.orbit(representation='trivial'))
        (1,2,3)
        >>> print(s.orbit(representation='sign'))
        - (1,2,3)

        '''
        if self.convention == 'Berger-Fresse':
            answer = self.zero()
            for k, v in self.items():
                seen = []
                for i in k:
                    if i not in seen:
                        seen.append(i)
                inverse = tuple(seen.index(i + 1) + 1 for i in range(len(seen)))
                inverse = SymmetricGroup_element(inverse)
                permutation = SymmetricModule_element({inverse: 1},
                                                      torsion=self.torsion)
                if representation == 'sign':
                    permutation = inverse.sign * permutation
                answer += permutation * self.create({k: v})

            return answer

        if self.convention == 'McClure-Smith':
            raise NotImplementedError

    def _reduce_rep(self):

        # remove non-surjections
        zeros = list()
        for k in self.keys():
            if set(k) != set(range(1, max(k) + 1)):
                zeros.append(k)
        for k in zeros:
            del self[k]

        # removes keys w/ equal consecutive values
        for k, v in self.items():
            for i in range(len(k) - 1):
                if k[i] == k[i + 1]:
                    self[k] = 0

        super()._reduce_rep()

    def table_arrangement(surj, only_dict=False):
        '''Returns the table arrangement of a surjection, as a tuple.
        If only_dict=True, it returns a dictionary dict such that
        dict[index] is the table arrangement of surj where
        surj[index] lies. surj must be a tuple or a list.
        '''

        def _final_indices(surj):
            '''Return the set of indices of elements in surj that are
               the last occurence of a value. surj must be a tuple or a list.
            '''
            finals = set()
            for k in range(1, max(surj) + 1):
                for index in reversed(range(len(surj))):
                    if surj[index] == k:
                        finals.add(index)
                        break
            return finals

        finals = _final_indices(surj)

        if only_dict:
            result = dict()
            row = 0
            for index, el in enumerate(surj):
                result[index] = row
                if index not in finals or index == len(surj) - 1:
                    row += 1
            return result

        result = tuple()
        new_row = tuple()
        for index, el in enumerate(surj):
            new_row += (el,)
            if index not in finals or index == len(surj) - 1:
                result += (new_row,)
                new_row = tuple()
        return result

    def last_index(value, sequence):
        '''returns the index of the last occurence of a value in a sequence'''

        for index, element in reversed(list(enumerate(sequence))):
            if value == element:
                return index
        raise ValueError(f"{value} is not an element of {sequence}")

    def degrees(surj, pieces):
        '''returns a tuple with the degree of the pieces in surj.
           surj must be a tuple or a list.

        '''
        row = Surjection_element.table_arrangement(surj, only_dict=True)
        result = tuple()
        remaining_surj = surj
        for piece in reversed(pieces):
            end = piece[-1]
            end_index = Surjection_element.last_index(end, remaining_surj)
            start_index = end_index - len(piece) + 1
            degree = row[end_index] - row[start_index]
            result = (degree,) + result
            remaining_surj = remaining_surj[0:start_index + 1]
        return result

    def _pcompose(self, other, k):
        '''partial composition self o_k other

        >>> u = Surjection_element({(1, 2, 1, 3): 1})
        >>> v = Surjection_element({(1, 2, 1): 1})
        >>> print(u._pcompose(v, 1))
        - (1,2,1,3,1,4) + (1,3,1,2,1,4) - (1,2,3,2,1,4)

        >>> m = Surjection_element({(1, 2): 1}, torsion=2)
        >>> br = Surjection_element({(1, 2, 1, 3, 1): 1}, torsion=2)
        >>> print(br.compose(m, 1))
        (1,3,1,2,4,2) + (1,2,3,2,4,2) + (1,3,1,4,1,2)

        # chain map check

        >>> u = Surjection_element({(1, 2, 1, 3): 1})
        >>> du = u.boundary()
        >>> v = Surjection_element({(1, 2, 1): 1})
        >>> dv = v.boundary()
        >>> du_v = du.compose(v, 1)
        >>> u_dv = u.compose(dv, 1)
        >>> uv = u._pcompose(v, 1)
        >>> duv = uv.boundary()
        >>> du_v - u_dv == duv
        True

        '''
        answer = self.zero()
        for k1, coeff1 in self.items():
            for k2, coeff2 in other.items():
                occurences = tuple(i for i, el in enumerate(k1) if el == k)

                k2_shifted = tuple(el + k - 1 for el in k2)
                k1_shifted = tuple(map(lambda i: i + other.arity - 1
                                       if i > k else i, k1))
                for k2_cut in decompositions(k2_shifted, len(occurences) - 1):
                    k3 = list(k1_shifted)
                    for x, idx in reversed(list(zip(k2_cut, occurences))):
                        k3[idx: idx + 1] = x

                    sign = 1
                    if self.arity != 2:  # Signs done by others, check
                        indices_to_cut_k1 = (0,) + occurences + (len(k1) - 1,)
                        k1_cut = tuple()
                        for i in range(len(indices_to_cut_k1) - 1):
                            left = indices_to_cut_k1[i]
                            right = indices_to_cut_k1[i + 1]
                            k1_cut += (k1_shifted[left:right + 1],)
                        deg1 = Surjection_element.degrees(k1_shifted, k1_cut)
                        deg2 = Surjection_element.degrees(k2_shifted, k2_cut)
                        sign_exp = 0
                        for index2 in range(len(k2_cut)):
                            deg_piece2 = deg2[index2]
                            left = occurences[index2]
                            right = len(k1_cut)
                            deg_pieces1 = sum(deg1[index1]
                                              for index1 in range(left, right))
                            sign_exp += (deg_piece2 * deg_pieces1) % 2
                        sign = (-1)**sign_exp

                    new_coeff = coeff1 * coeff2 * sign
                    answer += self.create({tuple(k3): new_coeff})
        return answer

    def compose(self, *others):
        '''general composition of Surjection_elements'''

        # partial composition
        if len(others) == 2 and isinstance(others[1], int):
            # unpaking and checking input
            other, k = others
            if self.torsion != other.torsion:
                raise TypeError('not the same torsion')
            return self._pcompose(other, k)
        # total composition
        else:
            if not len(others) == self.arity:
                raise TypeError('the number of arguments must be equal to '
                                + 'the arity of self')
            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer._pcompose(other, idx + 1)
            return answer

    def interval_cut(self, n):
        '''...'''
        def all_cuts(k, n):
            '''(k,n) -> 0 = n0 <= n1 <= ... <= nk <= n'''
            if k == 0:
                yield (0,)
            if k > 0:
                for cut in all_cuts(k - 1, n):
                    for i in range(cut[-1], n + 1):
                        yield cut + (i,)

        def constraints(surj):
            '''...'''
            for i in range(1, max(surj) + 1):
                preimage = (idx for idx, s in enumerate(surj) if s == i)
                for s, t in pairwise(preimage):
                    yield (s + 1, t)

        def good_cut(cut, const):
            '''...'''
            return all(cut[i] != cut[ipp] for i, ipp in const)

        def cut2multioperator(cut, n):
            '''...'''
            multisimplex = {i: tuple() for i in range(1, max(surj) + 1)}
            for i, pair in enumerate(pairwise(cut)):
                multisimplex[surj[i]] += tuple(range(pair[0], pair[1] + 1))

            std = set(range(n + 1))  # {0,1,...,n}
            return tuple(tuple(std.difference(set(spx)))
                         for spx in multisimplex.values())

        def cut_sign(cut, surj):

            class LeveledInterval:
                def __init__(self, start, end, level):
                    self.start = start
                    self.end = end
                    self.level = level
                    self.is_inner = True

                def length(self):
                    length = interval.end - interval.start
                    if self.is_inner:
                        length += 1
                    return length

            # transform cut to tuple of intervals
            intervals = tuple(LeveledInterval(cut[i], cut[i + 1], surj[i])
                              for i in range(len(surj)))

            # classify intervals: internal & final
            for level in range(1, max(surj) + 1):
                for interval in reversed(intervals):
                    if level == interval.level:
                        interval.is_inner = False
                        break

            # position sign = sum of n_i over internal intervals (n_{i-1}, n_i)
            position_sign_exp = sum(interval.end for interval in intervals
                                    if interval.is_inner)

            # permutation sign associated to ordering the surjection
            ordered = sorted(intervals, key=attrgetter('level'))
            perm_sign_exp = sum(ordered[i].length() * ordered[j].length()
                                for i in range(len(ordered))
                                for j in range(i, len(ordered))
                                if intervals.index(ordered[i])
                                > intervals.index(ordered[j]))

            return (-1)**((position_sign_exp + perm_sign_exp) % 2)

        answer = EilenbergZilber_element().copy_attrs_from(self)
        for surj, coeff in self.items():
            const = set(constraints(surj))
            k = len(surj) - 1
            good_cuts = (cut + (n,) for cut in all_cuts(k, n)
                         if good_cut(cut + (n,), const))
            for cut in good_cuts:
                sign = cut_sign(cut, surj)
                multiop = cut2multioperator(cut, n)
                answer += EilenbergZilber_element(
                    {multiop: sign * coeff}).copy_attrs_from(answer)
        return answer

    def _index_occurence(value, n, sequence):
        '''returns the index of the n-th occurence of value in sequence
        '''
        count = 0
        for index, element in enumerate(sequence):
            if element == value:
                count += 1
                if count == n:
                    return index

    def tau_vertex(vertex, surj):
        '''...'''
        indices = [Surjection_element._index_occurence(val + 1, x + 1, surj)
                   for val, x in enumerate(vertex)]
        indices.sort()
        return tuple(surj[i] for i in indices)

    def prism(surj):
        '''...'''
        r = max(surj)
        num_occurences = [surj.count(k) for k in range(1, r + 1)]

        # compute the base prism {0, ..., d1 - 1} x ... x {0, ..., dr - 1}
        summands = [range(num_occurences[k]) for k in range(r)]

        base_prism = product(*summands)
        # better to make summands a generator instead of using * ?

        base_prism = tuple(vertex for vertex in base_prism)
        return base_prism

    def tau(surj):
        '''for prismatic decomposition'''

        base_prism = Surjection_element.prism(surj)

        image_prism = tuple(Surjection_element.tau_vertex(vertex, surj)
                            for vertex in base_prism)
        return base_prism, image_prism

    def max_simplex(vect, r):
        '''maximal simplex determined by the sequence
           vect = (k0, ..., kd)
        '''

        result = tuple()
        vertex = tuple(0 for k in range(r))
        result = (vertex,)

        for k in vect:
            k = k - 1
            vertex = vertex[:k] + (vertex[k] + 1,) + vertex[k + 1:]
            result = result + (vertex,)

        return result

    def table_completion(self):
        '''...'''

        result = BarrattEccles_element(torsion=self.torsion)
        r = self.arity

        for surj, coeff in self.items():
            finals = Surjection_element._final_indices(surj)
            fund_vect = tuple(el for i, el in enumerate(surj)
                              if i not in finals)
            fund_simplex = Surjection_element.max_simplex(fund_vect, r)
            tau_fund_simplex = tuple(Surjection_element.tau_vertex(vx, surj)
                                     for vx in fund_simplex)

            # make a list with all the possibilities for indices k_i,
            # including eventual repetitions
            possibilities = []
            for k in range(1, r + 1):
                num_occurences = surj.count(k)
                possibilities += [k] * (num_occurences - 1)

            # iterate over all tuples (k_0, ..., k_d), where each k_i
            # has as many repetitions as it has in `possibilities`
            # Note: distinct_permutations need the more_itertools module.
            # We can replace it by ´set(permutations(possibilities))´,
            # but it will compute much more than needed.
            for vect in distinct_permutations(possibilities):
                simplex = Surjection_element.max_simplex(vect, r)

                # compute the image of the simplex
                tau_simplex = tuple(Surjection_element.tau_vertex(vx, surj)
                                    for vx in simplex)

                # compute the sign of the simplex
                sgn = 1

                for index, vertex in enumerate(tau_simplex):
                    vertex_f = tau_fund_simplex[index]

                    # find the permutation that takes vertex to vertex_f
                    permutation = {}
                    for idx, el in enumerate(vertex):
                        permutation[idx] = vertex_f.index(el)

                    # compute the sign of the permutation
                    sgn_perm = 1
                    for i in range(len(vertex)):
                        for j in range(i + 1, len(vertex)):
                            diff = permutation[j] - permutation[i]
                            sgn_perm *= diff // abs(diff)

                    sgn *= sgn_perm

                result += BarrattEccles_element({tau_simplex: coeff * sgn},
                                                torsion=self.torsion)
        return result

    def __call__(self, other):
        '''...'''
        assert self.torsion == other.torsion, "defined for the same ring"
        assert isinstance(self.degree, int) and self.arity, "defined for homogeneous surjections"

        if isinstance(other, CubicalEilenbergZilber_element):
            assert other.arity == 1, "defined for chains on a single cube"
            answer = other.zero()
            iterated_diagonal = other
            for _ in range(self.degree + self.arity - 1):
                iterated_diagonal = iterated_diagonal.coproduct()
            for k1, v1 in self.items():
                for k2, v2 in iterated_diagonal.items():
                    # sign
                    odds = [i for i, x in enumerate(k2) if x.count('e') % 2]
                    coeff = v1 * v2
                    for idx, i in enumerate(odds):
                        coeff *= (-1)**len([j for j in odds[idx + 1:] if k1[i] > k1[j]])
                    # elements
                    elements = []
                    for s in range(1, max(k1) + 1):
                        element = other.create({tuple(k2[i] for i, s_i in enumerate(k1)
                                                      if s_i == s): 1})
                        elements.append(element.product())
                    if all(elements):
                        for multipair in product(*(element.items() for element in elements)):
                            new_key = tuple(pair[0][0] for pair in multipair)
                            answer += answer.create({new_key: coeff})
            return answer


class Surjection():
    '''Class producing Surjection elements of special interest.'''

    @staticmethod
    def steenrod_product(arity, degree, torsion=None,
                         convention='Berger-Fresse'):
        '''Returns a surjection element representing the Steenrod
        product in the given arity and degree.

        Constructed recursively by mapping the minimal resolution W(r) 
        of Z[S_r] to Surj(r). We use the chain homotopy equivalence 
        of Surj(r) and Z defined using the chain contraction (i, p, s) 
        relating Surj(r-1) and Surj(r).

        Parameters
        ----------
        arity : int 
        Arity of the complex considered, Surj(arity).
        
        degree : int
        degree of the element considered Surj(arity)_degree.


        Examples
        --------

        # chain map checks:

        >>> t = SymmetricModule.transposition_element(6)
        >>> x = Surjection.steenrod_product(6, 3).boundary()
        >>> y = t * Surjection.steenrod_product(6, 2)
        >>> print(x == y)
        True

        >>> n = SymmetricModule.norm_element(5, torsion=7)
        >>> x = Surjection.steenrod_product(5, 6, torsion=7).boundary()
        >>> y = n * Surjection.steenrod_product(5, 5, torsion=7)
        >>> print(x == y)
        True

        >>> t = SymmetricModule.transposition_element(5)
        >>> x = Surjection.steenrod_product(5, 3, \
            convention='McClure-Smith').boundary()
        >>> y = t * Surjection.steenrod_product(5, 2, \
            convention='McClure-Smith')
        >>> print(x == y)
        True


        '''

        def i(surj, iterate=1):
            '''Inclusion of Surj(r) into Surj(r+1) by appending 1 at
            the start of basis elements and raising the value of all
            other entries by 1.'''

            if iterate == 1:
                answer = surj.zero()
                for k, v in surj.items():
                    answer += answer.create(
                        {(1,) + tuple(j + 1 for j in k): v})
                return answer
            if iterate > 1:
                return i(i(surj, iterate=iterate - 1))

        def p(surj, iterate=1):
            '''Projection of Surj(r) to Surj(r-1) by removing 1
            from a basis element with only one occurence of value 1
            and substracting 1 from all other entries.

            '''
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
            '''Chain homotopy from the identity to the composition pi.'''

            answer = surj.zero()
            for k, v in surj.items():
                answer += answer.create({(1,) + tuple(j for j in k): v})
            return answer

        def h(surj):
            '''Chain homotopy from i...i p..p to the identiy in Surj(r)
            realizing its contractibility to Surj(1).

            '''
            answer = s(surj)
            for r in range(1, arity - 1):
                answer += i(s(p(surj, r)), r)
            return answer

        operators = {
            0: SymmetricModule.norm_element(arity),
            1: SymmetricModule.transposition_element(arity)
        }

        def _psi(arity, degree, convention=convention):
            '''Recursive definition over the integers.'''

            if degree == 0:
                return Surjection_element({tuple(range(1, arity + 1)): 1},
                                          convention=convention)
            else:
                previous = _psi(arity, degree - 1, convention=convention)
                acted_on = operators[degree % 2] * previous
                answer = h(acted_on)
                return answer

        integral_answer = _psi(arity, degree, convention=convention)
        if torsion:
            integral_answer.set_torsion(torsion)
        return integral_answer

    @staticmethod
    def basis(arity, degree, complexity=None):
        ''' Returns the list of tuples forming the basis of the surjection
        operad in the given degree, arity and complexity

        >>> basis = sorted(Surjection.basis(3, 5, 3))
        >>> for s in basis:
        ...     print(s)
        (1, 2, 1, 3, 1, 3, 2, 3)
        (1, 3, 1, 2, 1, 2, 3, 2)
        (2, 1, 2, 3, 2, 3, 1, 3)
        (2, 3, 2, 1, 2, 1, 3, 1)
        (3, 1, 3, 2, 3, 2, 1, 2)
        (3, 2, 3, 1, 3, 1, 2, 1)
        '''

        if complexity is None:
            complexity = degree + 1
        a, d, c = arity, degree, complexity
        basis = []
        for s in product(range(1, a + 1), repeat=a + d):
            surj = Surjection_element({s: 1})
            if surj and surj.complexity <= c and surj.arity == a:
                basis.append(s)

        return basis


if __name__ == "__main__":
    import doctest
    doctest.testmod()
