def _cubical_action(self, other):
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
