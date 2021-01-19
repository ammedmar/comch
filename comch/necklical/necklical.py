from ..free_module import FreeModuleElement
from ..simplicial import Simplex, SimplicialElement
from functools import reduce

# a method in the class Cube
def necklace(self):
    """..."""
    ones = [idx for idx, i in enumerate(self) if i == 1]
    twos = [idx for idx, i in enumerate(self) if i == 2]
    aux = [0]
    answer = []
    for idx in range(len(self) + 2):
        if idx in twos:
            aux.append(idx + 1)
        if idx in ones:
            aux.append(idx + 1)
            answer.append(Simplex(aux))
            aux = [idx + 1]
    if aux:
        aux.append(len(self) + 1)
        answer.append(Simplex(aux))
    return Necklace(answer)

class Necklace(tuple):
    """Tuple of simplices."""

    def __init__(self, simplices):
        answer = tuple(Simplex(spx) for spx in simplices)

    def one_reduced(self):
        return Necklace(filter(lambda spx: spx.dimension != 1, self))


class NecklicalElement(FreeModuleElement):
    """..."""

    def __init__(self, data=None, torsion=None):
        if data:
            new_data = {}
            for k, v in data.items():
                new_k = tuple(Necklace(necklace) for necklace in k)
                new_data[new_k] = v
            data = new_data

        super(NecklicalElement, self).__init__(
            data=data, torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def one_reduced(self):
        answer = self.zero()
        for k, v in self.items():
            answer += answer.create({tuple(n.one_reduced() for n in k): v})
        return answer

    def to_simplicial(self, one_reduced=True):
        answer = SimplicialElement(torsion=self.torsion)
        for k, v in self.items():
            new_k = tuple(reduce(lambda x, y: x + y, k))
            if one_reduced:
                new_k = tuple(spx for spx in new_k if len(spx) > 2)
            answer += answer.create({new_k: v})
        return answer
