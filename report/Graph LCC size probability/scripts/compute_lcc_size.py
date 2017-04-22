#!/usr/bin/python3

from time import time
from math import factorial  # effective enough, might implement tree factorial algorithm
                            # if computation time gets too large

from PermWeightedSum import PermWeightedSum
from Combinatorics import Combinatorics
from ConnectedCounter import ConnectedCounter

def X(n):
    '''return n choose 2'''
    return n*(n-1)//2

def Pkalpha(n, k, alpha):
    return (factorial(n)) // (factorial(alpha) * (factorial(k))**alpha * factorial(n - k*alpha))

def get_permutations_that_sum_to(value, nb_elements, beg=0, end=-1):
    if end < 0:
        end=value
    if end < beg or nb_elements * beg > value or nb_elements * end < value:
        raise StopIteration
    perm = [beg] * nb_elements
    # Can be optimized...
    while True:
        if sum(perm) == value:
            yield perm
        perm[-1] += 1
        idx = nb_elements-1
        while idx >= 1 and perm[idx] == end+1:
            perm[idx] = beg
            idx -= 1
            perm[idx] += 1
        if perm[0] > end:
            break

class GSS:  # Graph Set Size
    # TODO: implement a cache

    @staticmethod
    def connected(n, m):
        return int(ConnectedCounter.count(n, m))

    @staticmethod
    def lcc_k(n, k):
        ret = 0
        for m in range(k-1, X(n)+1):
            ret += GSS.lcc_km(n, k, m)
        return ret

    @staticmethod
    def lcc_km(n, k, m):
        if k == 0:
            if n == m == 0:
                return 1
            else:
                return 0
        ret = 0
        for alpha in range(1, n//k + 1):
            ret += GSS.lcc_kmalpha(n, k, m, alpha)
        return ret

    @staticmethod
    def lcc_kmalpha(n, k, m, alpha):
        if k == 0 and (m != 0 or n != 0):
            return 0
        if alpha > n//k or m < k-1 or m > n*(k-1)/2 or k > n:
            return 0
        if m == 0:
            return int(
                (k == alpha == 1)
                    or
                (k == n == m)
            )
        if m == 1:
            return X(n) if k == 2 and alpha == 1 else 0
        if k == n:  # and alpha == 1
            return GSS.connected(n, m)
        ret = 0
        for Sigma in range(alpha*(k-1), min(m, alpha*X(k))+1):
            for perm in get_permutations_that_sum_to(Sigma, alpha, beg=k-1, end=X(k)):
                prod = 1
                for j in range(alpha):
                    prod *= GSS.connected(k, perm[j])
                tmp = 0
                for p in range(k):
                    tmp += GSS.lcc_km(n-k*alpha, p, m-Sigma)
                ret += tmp*prod
        return ret * Pkalpha(n, k, alpha)

##### tests

def test_connected():
    for n in range(2, 7):
        for k in range(n-1, X(n)+1):
            print('(n, k) == {} yields {}'.format((n, k), GSS.connected(n, k)))
        print('')

def test_Pkalpha():
    print(Pkalpha(6, 3, 2))
    print(Pkalpha(7, 2, 3))

def test_permutations():
    for perm in get_permutations_that_sum_to(6, 3, 1, 3):
        print(perm)
    for perm in get_permutations_that_sum_to(10, 4):
        print(perm)

def test_lcc_km():
    for m in range(3):
        print('(n, k, m) == {}   |--> {}'.format((5, 2, m), GSS.lcc_km(5, 2, m)))
        print('expected: {}\n\n'.format(Pkalpha(5, 2, m) if m > 0 else 0))
    total = 0
    N = 5
    m_cache = [0] * (X(N)+1)
    for k in range(1, N+1):
        for m in range(k-1, X(N)+1):
            tmp = GSS.lcc_km(N, k, m)
            m_cache[m] += tmp
            total += tmp
            print('\t(n, k, m) == {}  --> {}'.format((N, k, m), tmp))
        print('')
    print('total is: {}  --- expected: {}'.format(total, 2**X(N)))
    for (m, v) in enumerate(m_cache):
        print('for m == {}, there are {} graphs --- expected: {}'.format(m, v, Combinatorics.binom(X(N), m)))

def test_sum_lcc_k():
    N = 5
    total = 0
    for k in range(1, N+1):
        total += GSS.lcc_k(N, k)
    print('total: {}  ---  expected: {}'.format(total, 2**X(N)))

    for N in range(3, 30):
        print(N)
        print(sum([GSS.lcc_k(N, k) for k in range(1, N+1)]), 2**X(N),sep='\n', end='\n\n')

def test_sum_connected():
    for N in range(3, 30):
        print('{}\n{}\n'.format(N, sum([GSS.connected(N, k) for k in range(N-1, X(N)+1)])))

if __name__ == '__main__':
    #test_connected()
    #test_Pkalpha()
    #test_permutations()
    test_lcc_km()
    #print('\n\n\t----------------\n\n')
    test_sum_lcc_k()
    #test_sum_connected()
