#!/usr/bin/env python

# this file is part of the github repository: https://github.com/nwhoppe/rosalind
# author: nwhoppe
# created: 9/3/18

import argparse


def mortal_fibonacci_rabbits(n, m):
    """
    calculate number of rabbit pairs at month n assuming one month to reach reproductive age and life span of m months.

    Args:
        n: 0 <= integer <= 100
        m: 0 <= integer <= 20

    Returns:
         The number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.

    Raises:
        ValueError: n or m are out of bounds.
    """

    if n < 0 or n > 100:
        raise ValueError("n must be in range 0 <= n <= 100\n")
    if m < 0 or m > 20:
        raise ValueError("m must be in range 0 <= m <= 20\n")
    if n == 0 or m == 1:
        return 0
    if n == 1 or n == 2 or m == 2:
        return 1

    rabbit_pairs = [0, 1]
    for i in xrange(2, n + 1):
        if i <= m:
            # before rabbits die, the number of pairs is equal to the standard fibonacci number
            rabbit_pairs.append(rabbit_pairs[i - 1] + rabbit_pairs[i - 2])
        elif i == m + 1:
            # initial rabbit pair death
            rabbit_pairs.append(rabbit_pairs[i - 1] + rabbit_pairs[i - 2] - 1)
        else:
            # rabbits pairs are now F_n-1 + F_n-12 - F_n-(m+1)
            rabbit_pairs.append(rabbit_pairs[i - 1] + rabbit_pairs[i - 2] - rabbit_pairs[i - (m + 1)])
    return rabbit_pairs[n]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to calculate rabbit-pair population at month n assuming 
    rabbits take one month to reach reproductive age and live m months before dying""")
    required = parser.add_argument_group('required')
    required.add_argument("-n", type=int, help="month number to calculate rabbit population for")
    required.add_argument("-m", type=int, help="rabbit live span in months")
    args = parser.parse_args()
    print mortal_fibonacci_rabbits(args.n, args.m)
