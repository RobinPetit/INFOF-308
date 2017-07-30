from Combinatorics import Combinatorics

class ConnectedCounter:
    '''SEE MSE #689526'''
    cache = [[-1 for m in range(n*(n-1)//2-(n-1)+1)] for n in range(100+1)]

    @staticmethod
    def count(*kwargs):
        return ConnectedCounter.qq(*kwargs)

    @staticmethod
    def qq(n, k):
        if k < n-1 or k > n*(n-1)//2:
            return 0
        else:
            if ConnectedCounter.cache[n][k-(n-1)] >= 0:
                return ConnectedCounter.cache[n][k-(n-1)]
            res = Combinatorics.binom(n*(n-1)//2, k)
            for m in range(n-1):
                res1 = 0
                for p in range(max(0, k-m*(m+1)//2), k-m+1):
                    res1 += Combinatorics.binom((n-1-m)*(n-2-m)//2, p) * ConnectedCounter.qq(m+1, k-p)
                res -= Combinatorics.binom(n-1, m) * res1
            ConnectedCounter.cache[n][k-(n-1)] = res
            return res
