from numpy import *  #pi, exp, sqrt, etc.
from scipy.special import erf, erfc


def gaussAssymetry(x, mu, sig, alpha):
    model = 1+erf((alpha*(x-mu))/(sig*sqrt(2)))
    return model/2

def lorentzAssymetry(x, mu, sig, alpha):
    model = 1+arctan(alpha*(x-mu)/sig)/(pi/2)
    return model/2