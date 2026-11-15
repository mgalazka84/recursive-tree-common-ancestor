"""Enumerate recursive trees and check the full joint point probabilities."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, factorial


def record_pmf(first, last):
    """Exact distribution of independent Bernoulli(1/i) sums; empty sum is zero."""
    numerators = [1]
    denominator = 1
    for i in range(first, last + 1):
        numerators = [(i - 1) * same + prev
                      for same, prev in zip(numerators + [0], [0] + numerators)]
        denominator *= i
    return [Fraction(x, denominator) for x in numerators]


def tree_counts(n):
    counts = Counter()
    for parents in product(*(range(child) for child in range(1, n))):
        depths = [0]
        ancestor_masks = [1]
        for child, parent in enumerate(parents, 1):
            depths.append(depths[parent] + 1)
            ancestor_masks.append(ancestor_masks[parent] | (1 << child))
        for u in range(n):
            for v in range(u + 1, n):
                ancestor = (ancestor_masks[u] & ancestor_masks[v]).bit_length() - 1
                h = depths[ancestor]
                left, right = depths[u] - h, depths[v] - h
                k = left + right
                counts[ancestor + 1, k, left, h] += 1
                counts[ancestor + 1, k, right, h] += 1
    return counts


def verify(n):
    total = factorial(n - 1) * n * (n - 1)
    counts = tree_counts(n)
    assert sum(counts.values()) == total
    if n == 7:
        distance_two = sum(count for (_, distance, _, _), count in counts.items()
                           if distance == 2)
        through_root = sum(count for (ancestor, distance, _, _), count in counts.items()
                           if distance == 2 and ancestor == 1)
        root_probability = Fraction(through_root, distance_two)
        assert (through_root, distance_two) == (6496, 10224)
        assert root_probability == Fraction(406, 639)
        print(f'  Distance-two example: {through_root:,}/{distance_two:,} = {root_probability}')
    for m in range(1, n + 1):
        suffix = record_pmf(m + 1, n)
        height = record_pmf(1, m - 1)
        for k in range(1, n):
            distance_and_label = (Fraction(2**k, (n - 1) * m) *
                                  (suffix[k] if k < len(suffix) else 0))
            for arm in range(k + 1):
                split = Fraction(comb(k, arm), 2**k)
                for h, hg in enumerate(height):
                    observed = counts.pop((m, k, arm, h), 0)
                    expected = total * distance_and_label * split * hg
                    assert observed == expected, (n, m, k, arm, h, observed, expected)
    assert not counts
    print(f'n={n}: {factorial(n - 1):,} trees, {total:,} ordered pairs, exact agreement')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--max-n', type=int, default=9)
    args = parser.parse_args()
    if args.max_n < 2:
        parser.error('--max-n must be at least 2')
    for n in range(2, args.max_n + 1):
        verify(n)
