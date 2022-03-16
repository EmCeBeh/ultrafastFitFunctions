from numpy import *  #pi, exp, sqrt, etc.
from scipy.special import erf, erfc

from .gauss import normGauss
from .assymetries import gaussAssymetry

def absorptionPeak(x, mu, sig, asym_plateau, asym_peak,  A_plateau, A_peak):

    plateau = A_plateau*gaussAssymetry(x, mu, 10, asym_plateau)
    peak    = A_peak*normGauss(x, mu, sig, 1, 0)*gaussAssymetry(x, mu, sig, asym_peak)
    
    return plateau + peak