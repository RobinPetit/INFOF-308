import math
import numpy as np

BINOM_CACHE_SIZE = 100
FACTORIAL_CACHE_SIZE = 1000

class Combinatorics:
    binom_cache = [[0 for k in range(n)] for n in range(BINOM_CACHE_SIZE+1)]
    factorial_cache = [1 for n in range(FACTORIAL_CACHE_SIZE+1)]
    factorial_last_cache = 0

    @staticmethod
    def binom(n, k):
        if k > n or k < 0:
            return 0
        if k == n or k == 0:
            return 1
        if k == 1 or k == n-1:
            return n
        if n <= BINOM_CACHE_SIZE and Combinatorics.binom_cache[n-1][k-1] > 0:
            return Combinatorics.binom_cache[n-1][k-1]
        ret = Combinatorics.factorial(n) // Combinatorics.factorial(k) // Combinatorics.factorial(n-k)
        if n <= BINOM_CACHE_SIZE:
            Combinatorics.binom_cache[n-1][k-1] = ret
        return ret

    @staticmethod
    def multinom(n, alpha):
        # assert sum(alpha) == n
        ret = Combinatorics.factorial(n)
        for el in alpha:
            ret //= Combinatorics.factorial(el)
        return ret

    @staticmethod
    def factorial(n):
        if n > FACTORIAL_CACHE_SIZE:
            return math.factorial(n)
        if n <= Combinatorics.factorial_last_cache:
            return Combinatorics.factorial_cache[n]
        ret = Combinatorics.factorial_cache[Combinatorics.factorial_last_cache]
        for i in range(Combinatorics.factorial_last_cache+1, n+1):
            ret *= i
            Combinatorics.factorial_cache[i] = ret
        Combinatorics.factorial_last_cache = n
        return ret

def test_factorial():
    from timeit import repeat
    NB_TESTS = 5000
    print(repeat("for i in range(1000+1):\n\tCombinatorics.factorial(i)", "from __main__ import Combinatorics", number=NB_TESTS))
    print(repeat("for i in range(1000+1):\n\tmath.factorial(i)", "import math", number=NB_TESTS))
    print(repeat("for i in range(1001, -1, -1):\n\tCombinatorics.factorial(i)", "from __main__ import Combinatorics", number=NB_TESTS))
    print(repeat("for i in range(1001, -1, -1):\n\tmath.factorial(i)", "import math", number=NB_TESTS))
    # gives following result:
    # [1.2725156589999642, 1.2541886670001077, 1.2880517129997315]
    # [63.86013914499972, 64.1316873830001, 66.46548106799992]
    # [1.450911523999821, 1.4265047030003188, 1.4306807879997905]
    # [64.40737520200037, 64.23296190400015, 68.07413684100038]


if __name__ == '__main__':
    test_factorial()
