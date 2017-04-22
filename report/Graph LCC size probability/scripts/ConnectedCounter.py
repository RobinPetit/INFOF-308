from Combinatorics import Combinatorics
from PermWeightedSum import PermWeightedSum

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
            if ConnectedCounter.is_in_cache(n, k):
                return ConnectedCounter.get_from_cache(n, k)
            res = Combinatorics.binom(n*(n-1)//2, k)
            for m in range(n-1):
                res1 = 0
                for p in range(max(0, k-m*(m+1)//2), k-m+1):
                    res1 += Combinatorics.binom((n-1-m)*(n-2-m)//2, p) * ConnectedCounter.qq(m+1, k-p)
                res -= Combinatorics.binom(n-1, m) * res1
            ConnectedCounter.set_in_cache(n, k, res)
            return res

    @staticmethod
    def set_in_cache(n, k, value):
        '''insert the given value in cache'''
        ConnectedCounter.cache[n][k-(n-1)] = value

    @staticmethod
    def get_from_cache(n, k):
        '''return the cached value for given graph dimensions'''
        return ConnectedCounter.cache[n][k-(n-1)]
        
    @staticmethod
    def is_in_cache(n, k):
        '''return True if value is in cache'''
        return ConnectedCounter.get_from_cache(n, k) > 0

