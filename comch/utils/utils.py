from itertools import chain, accumulate, permutations, tee
from operator import add


def partitions(n, k, smallest_value=1, largest_value=None, ordered=False):
    """n is the integer to partition and k is the length of partitions.
    It returns all k tuples of integers greater or equal to smallest_value
    and less than or equal to largest_value that add up to n.
    If ordered == True it returns all tuples if False it returns those
    in non-decreassing order
    """

    if largest_value is None:
        largest_value = n

    def unordered_partitions(n, k, r=smallest_value, m=largest_value):
        if k == 1:
            if r <= n <= m:
                yield (n,)
            return
        for i in range(r, m + 1):
            for result in unordered_partitions(n - i, k - 1, i, m):
                yield (i,) + result

    if ordered:
        return chain.from_iterable(set(permutations(p)) for p
                                   in unordered_partitions(n, k))
    if not ordered:
        return unordered_partitions(n, k)


def decompositions(x, k):
    """Returns a set with all tuples
    (x[:i_1+1], x[i_1: i_2+1], ..., x[i_{k-1}: i_k+1], x[i_k:])
    with i_1 <= i_2 <= ... <= i_k

    >>> sorted(decompositions((0, 1), 2))
    [((0,), (0,), (0, 1)), ((0,), (0, 1), (1,)), ((0, 1), (1,), (1,))]

    """
    if isinstance(x, int):
        for part in partitions(x + k, k + 1, ordered=True):
            augmented_desuspension = (0,) + tuple(i - 1 for i in part)
            indices = tuple(accumulate(augmented_desuspension,
                                       add))
            pairs = tuple((i, j + 1) for i, j in pairwise(indices))
            yield pairs
    else:
        for pairs in decompositions(len(x), k):
            yield tuple(x[pair[0]: pair[1]] for pair in pairs)


def distinct_permutations(iterable, r=None):
    """Yield successive distinct permutations of the elements in *iterable*.

    >>> sorted(distinct_permutations([1, 0, 1]))
    [(0, 1, 1), (1, 0, 1), (1, 1, 0)]

    Equivalent to ``set(permutations(iterable))``, except duplicates are not
    generated and thrown away. For larger input sequences this is much more
    efficient.

    Duplicate permutations arise when there are duplicated elements in the
    input iterable. The number of items returned is
    `n! / (x_1! * x_2! * ... * x_n!)`, where `n` is the total number of
    items input, and each `x_i` is the count of a distinct item in the input
    sequence.

    If *r* is given, only the *r*-length permutations are yielded.

    >>> sorted(distinct_permutations([1, 0, 1], r=3))
    [(0, 1, 1), (1, 0, 1), (1, 1, 0)]

    """
    # Algorithm: https://w.wiki/Qai
    def _full(A):
        while True:
            # Yield the permutation we have
            yield tuple(A)

            # Find the largest index i such that A[i] < A[i + 1]
            for i in range(size - 2, -1, -1):
                if A[i] < A[i + 1]:
                    break
            #  If no such index exists, this permutation is the last one
            else:
                return

            # Find the largest index j greater than j such that A[i] < A[j]
            for j in range(size - 1, i, -1):
                if A[i] < A[j]:
                    break

            # Swap the value of A[i] with that of A[j], then reverse the
            # sequence from A[i + 1] to form the new permutation
            A[i], A[j] = A[j], A[i]
            A[i + 1:] = A[: i - size: -1]  # A[i + 1:][::-1]

    # Algorithm: modified from the above
    def _partial(A, r):
        # Split A into the first r items and the last r items
        head, tail = A[:r], A[r:]
        right_head_indexes = range(r - 1, -1, -1)
        left_tail_indexes = range(len(tail))

        while True:
            # Yield the permutation we have
            yield tuple(head)

            # Starting from the right, find the first index of the head with
            # value smaller than the maximum value of the tail - call it i.
            pivot = tail[-1]
            for i in right_head_indexes:
                if head[i] < pivot:
                    break
                pivot = head[i]
            else:
                return

            # Starting from the left, find the first value of the tail
            # with a value greater than head[i] and swap.
            for j in left_tail_indexes:
                if tail[j] > head[i]:
                    head[i], tail[j] = tail[j], head[i]
                    break
            # If we didn't find one, start from the right and find the first
            # index of the head with a value greater than head[i] and swap.
            else:
                for j in right_head_indexes:
                    if head[j] > head[i]:
                        head[i], head[j] = head[j], head[i]
                        break

            # Reverse head[i + 1:] and swap it with tail[:r - (i + 1)]
            tail += head[: i - r: -1]  # head[i + 1:][::-1]
            i += 1
            head[i:], tail[:] = tail[: r - i], tail[r - i:]

    items = sorted(iterable)

    size = len(items)
    if r is None:
        r = size

    if 0 < r <= size:
        return _full(items) if (r == size) else _partial(items, r)

    return iter(() if r else ((),))


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)
