import numpy as np

from math import pow

def linear_regression(x, y):
    """
    returns the m and p coefficients of y = mx + p
    """
    x = np.array(x)
    y = np.array(y)
    A = np.vstack([x, np.ones(len(x))]).T
    m, p = np.linalg.lstsq(A, y)[0]
    return m, p
