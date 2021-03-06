#!/usr/bin/python3

from time import time
from Combinatorics import Combinatorics
from ConnectedCounter import ConnectedCounter

import numpy as np

# Avoid computing factorials
def X(n):
    '''return n choose 2'''
    return n*(n-1)//2

def Pkalpha(n, k, alpha):
    '''return the number of sets {W_1, ..., W_alpha} such that all the
    W_i's are taken disjointed in a set of size n, and that all the W_i's
    are of size k'''
    return Combinatorics.factorial(n) // \
            (Combinatorics.factorial(alpha) \
            * Combinatorics.factorial(k)**alpha \
            * Combinatorics.factorial(n - k*alpha))

# Remark: by only generating such vectors st v_{i-1} <= v_i would
# allow to generate less of these, letting user get all of the permutations
def get_permutations_that_sum_to(value, nb_elements, beg=0, end=-1):
    '''generating function that yields all integer vectors of given size that
    sum to a given value'''
    if end < 0:
        end=value
    if end < beg or nb_elements * beg > value or nb_elements * end < value:
        raise StopIteration
    perm = [beg] * nb_elements
    s = beg*nb_elements  # keep track of the current value of the sum
    # Loop can be optimized...
    while True:
        if s > value:
            s += end + 1 - perm[-1]
            perm[-1] = end+1
        elif s == value:
            yield perm
            perm[-1] += 1
            s += 1
        else:
            perm[-1] += value-s
            s = value
        idx = nb_elements-1
        # rearrange vector so that each element is in {beg, ..., end}
        while idx >= 1 and perm[idx] > end:
            perm[idx-1] += 1
            s += 1 + beg - perm[idx]
            perm[idx] = beg
            idx -= 1
        if perm[0] > end:
            break

LCC_CACHE_SIZE = 250
class LccCounter:
    '''
        Static class made for counting graphs having given largest connected component
    '''
    # allocating cache
    print('creating cache')
    cache = [[[0 for m in range(k-1, n*(n-1)//2+1)] for k in range(n+1)] for n in range(LCC_CACHE_SIZE+1)]
    print('cache created')

    @staticmethod
    def connected(n, m):
        '''return the number of connected graphs having n vertices and m edges'''
        return int(ConnectedCounter.count(n, m))

    @staticmethod
    def lcc_k(n, k):
        '''return the amount of graphs having n vertices and a largest connected
        component of size k'''
        ret = 0
        for m in range(k-1, X(n)+1):
            ret += LccCounter.lcc_km(n, k, m)
        return ret

    @staticmethod
    def lcc_km(n, k, m):
        '''return the number of graphs having n vertices, m edges, and a
        largest connected component of size k'''
        if k == 0:
            if n == m == 0:
                return 1
            else:
                return 0
        if m < k-1 or m > n*(n-1)//2 or k > n:
            return 0
        if n <= LCC_CACHE_SIZE:
            if LccCounter.cache[n][k][m-k+1] > 0:
                return LccCounter.cache[n][k][m-k+1]
        ret = sum([LccCounter.lcc_kmalpha(n, k, m, alpha) for alpha in range(1, n//k+1)])
        if n <= LCC_CACHE_SIZE:
            LccCounter.cache[n][k][m-k+1] = ret
        return ret

    @staticmethod
    def lcc_kmalpha(n, k, m, alpha):
        '''return the number of graphs having n vertices, m edges, a largest
        connected component of size k, and alpha connected components of size k'''
        if k == 0 and (m != 0 or n != 0):
            return 0
        if alpha > n//k: #or m < k-1 or m > n*(k-1)/2 or k > n:
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
            return LccCounter.connected(n, m)
        ret = 0
        for Sigma in range(alpha*(k-1), min(m, alpha*X(k))+1):
            for perm in get_permutations_that_sum_to(Sigma, alpha, beg=k-1, end=X(k)):
                prod = 1
                for j in range(alpha):
                    prod *= LccCounter.connected(k, perm[j])
                tmp = 0
                for p in range(k):
                    tmp += LccCounter.lcc_km(n-k*alpha, p, m-Sigma)
                ret += tmp*prod
        return ret * Pkalpha(n, k, alpha)

    @staticmethod
    def get_distribution_lcc_m(n, m):
        '''return a vector [p_k]_k such that p_k is the probability that |LCC(\Gamma)| = k
        for \Gamma a random graph having n vertices and m edges'''
        ret = [LccCounter.lcc_km(n, k, m) for k in range(n+1)]
        s = sum(ret)
        for i in range(n+1):
            ret[i] /= s
        return ret

    @staticmethod
    def get_distribution_lcc(n):
        '''return a vector [p_k]_k such that p_k is the probability that |LCC(\Gamma)| = k
        for \Gamma a random graph having n vertices'''
        ret = [LccCounter.lcc_k(n, k) for k in range(n+1)]
        s = 1 << X(n)  # 2**X(n)
        for k in range(n+1):
            ret[k] /= s
        return ret

    @staticmethod
    def get_expectation_lcc(n):
        '''return the expectation of the LCC of a random graph having n vertices'''
        return LccCounter.weighted_sum(LccCounter.get_distribution_lcc(n))

    @staticmethod
    def get_expectation_lcc_m(n, m):
        '''return the expectation of the LCC of a random graph having n
        vertices and m edges'''
        return LccCounter.weighted_sum(LccCounter.get_distribution_lcc_m(n, m))

    @staticmethod
    def get_std_lcc(n):
        '''return the standard deviation of the LCC of a random graph having n
        vertices'''
        expectation = LccCounter.get_expectation_lcc(n)
        return LccCounter.weighted_sum([(p_k - expectation)**2 for p_k in LccCounter.get_distribution_lcc(n)])

    @staticmethod
    def get_std_lcc_m(n, m):
        '''return the standard deviation of the LCC of a random graph having n
        vertices and m edges'''
        expectation = LccCounter.get_expectation_lcc(n)
        return LccCounter.weighted_sum([(p_k - expectation)**2 for p_k in LccCounter.get_distribution_lcc_m(n, m)])

    @staticmethod
    def weighted_sum(vect):
        return sum([idx*value for idx, value in enumerate(vect)])


##### tests

def test_connected():
    for n in range(2, 7):
        for k in range(n-1, X(n)+1):
            print('(n, k) == {} yields {}'.format((n, k), LccCounter.connected(n, k)))
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
        print('(n, k, m) == {}   |--> {}'.format((5, 2, m), LccCounter.lcc_km(5, 2, m)))
        print('expected: {}\n\n'.format(Pkalpha(5, 2, m) if m > 0 else 0))
    total = 0
    N = 5
    m_cache = [0] * (X(N)+1)
    for k in range(1, N+1):
        for m in range(k-1, X(N)+1):
            tmp = LccCounter.lcc_km(N, k, m)
            m_cache[m] += tmp
            total += tmp
            print('\t(n, k, m) == {}  --> {}'.format((N, k, m), tmp))
        print('')
    print('total is: {}  --- expected: {}'.format(total, 2**X(N)))
    for (m, v) in enumerate(m_cache):
        print('for m == {}, there are {} graphs --- expected: {}'.format(m, v, Combinatorics.binom(X(N), m)))
    print('sum: {}'.format(sum(m_cache)))

def test_sum_lcc_k():
    N = 5
    total = 0
    for k in range(1, N+1):
        total += LccCounter.lcc_k(N, k)
    print('total: {}  ---  expected: {}'.format(total, 2**X(N)))

    a = time()
    for N in range(3, 17):
        #print(N)
        print(sum([LccCounter.lcc_k(N, k) for k in range(1, N+1)]), 2**X(N),sep='\n', end='\n\n')
    b = time()
    print('with caching, time == {} ms'.format(int((b-a)*1000)))

def test_sum_connected():
    for N in range(3, 30):
        print('{}\n{}\n'.format(N, sum([LccCounter.connected(N, k) for k in range(N-1, X(N)+1)])))

########## main

if __name__ == '__main__':
    #test_connected()
    #test_Pkalpha()
    #test_permutations()
    test_lcc_km()
    #print('\n\n\t----------------\n\n')
    #test_sum_lcc_k()
    #test_sum_connected()
    for m in range(66+1):
        print('expected LCC size in \\Gamma_{}(12, .): {}'.format(m, LccCounter.get_expectation_lcc_m(12, m)))
    #distribution = get_distribution_lcc(12)
    #for k in range(len(distribution)):
    #    print('probability that |LCC| of a graph in \\Gamma(12, .) is {} == {} %'.format(k, 100*distribution[k]))
    #print('(summing to {})'.format(sum(distribution)))
