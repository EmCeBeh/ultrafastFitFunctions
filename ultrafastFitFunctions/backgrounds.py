from numpy import pi, exp, sqrt, cos
from scipy.special import erf, erfc


def cBg(x, mu, A, c):
    '''
    Charge background for reflectivity measurements
    In the X-ray range 'x' is equivalent to q / hkl values and given in e.g. invers angstroms
    '''
    
    model = A*(x-mu)**(-4) + c
    return model


def benOcko1(x, qc, a ,b, c, d, sig):
    '''
    x    
    qc   where the background starts
    a    overall scaling
    b    fringes size
    c    fringes frequency
    d    offset
    sig  sigma of the exponential decay
    '''
    qPrime = sqrt(x**2 - qc**2)
    Rf = abs((x-qPrime)/(x+qPrime))**2
    R = Rf*a*(1+b*cos(x*c))*exp(-(x*sig**2)) + d
    return R


def benOcko2(x, qc, a ,b, c, d, sig):

    qPrime = sqrt(x**2 - qc**2)
    
    Rf = abs((x-qPrime)/(x+qPrime))**2
    R = Rf*(a**2 + b**2 + 2*a*b*cos(x*d))*exp(-x**2*sig**2)

    return R

#def reflectivityPseudoVoigt(x, mu1, mu2, sig1, sig2, A1, A2, c1, alpha2):
#        model1 = A*(x-mu1)**(-4) + c1
#        model2 = B*(alpha2 * 1/(1 + ((x-mu2)/sig2)**2) + (1-alpha2)*exp(-log(2)*((x-mu2)/sig2)**2))
#        return model1 + model2
# Same as adding the two above...


def oscillation(x, mu, A, p):
    term1 = -A*cos(1/p*2*pi*(x-mu)) + A
    return term1*step(x, mu)

def oscillationStep(x, mu, phase, A, p, Aexp, tau, offset):
    term1 = -A*cos(1/p*2*pi*(x-mu)) + A
    term2 = Aexp*exp(-(x-mu)/tau) - Aexp
    return (term1+term2)*step(x, mu)