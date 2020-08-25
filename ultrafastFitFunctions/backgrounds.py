from numpy import *  #pi, exp, sqrt, etc.
from scipy.special import erf, erfc




def cBg(x, mu, A, c):
    '''
    Charge background for reflectivity measurements
    In the X-ray range 'x' is equivalent to q / hkl values and given in e.g. invers angstroms
    '''
    
    model = A*(x-mu)**(-4) + c
    return model

def benOcko(x, qc,a ,b, c, d, sig):

    qPrime = sqrt(x**2 - qc**2)
    
    Rf = abs((x-qPrime)/(x+qPrime))**2
    
#     R = Rf*(a**2 + b**2 + 2*a*b*cos(x*d))*exp(-x**2*sigma**2)
    R = Rf*a*(1+b*cos(x*d))*exp(-(x*sig**2)) + c

    return R

#def reflectivityPseudoVoigt(x, mu1, mu2, sig1, sig2, A1, A2, c1, alpha2):
#        model1 = A*(x-mu1)**(-4) + c1
#        model2 = B*(alpha2 * 1/(1 + ((x-mu2)/sig2)**2) + (1-alpha2)*exp(-log(2)*((x-mu2)/sig2)**2))
#        return model1 + model2
# Same as adding the two above...