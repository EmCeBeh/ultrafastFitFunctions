from numpy import *  # pi, exp, sqrt, etc.
from scipy.special import erf, erfc


def gauss(x, sig):
    r"""
    What does the function do?

    Parameters
    ----------
    x : array_like
        `x` represents the x-coordinates of a set of datapoints.
    a : int
        `a` represents the variable, which does... .

    Returns
    -------
    y : ndarray
        An array of the same shape as x.

    """
    return 1 / sqrt(2 * pi) / sig * exp(-(x**2) / (2 * sig**2))


def gaussHeight(x, sig, A, c):
    return A * exp(-(x**2) / (2 * sigma**2)) + c


def normGauss(x, mu, sig, A, c):
    model = A * exp(-(((x - mu) / sig) ** 2) / 2) + c
    return model


def gaussAssymetry(x, mu, sig, alpha):
    """
    To be multiplied by gaussian above to provide assymetry
    """
    model = 1 + erf((alpha * (x - mu)) / (sig * sqrt(2)))
    return model
