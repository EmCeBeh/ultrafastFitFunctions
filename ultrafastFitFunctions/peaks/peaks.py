from numpy import *  #pi, exp, sqrt, etc.
from scipy.special import erf, erfc

def lorentz(x, mu, sig, A, c):
    model = A/(1 + ((x-mu)/sig)**2)
    return model


def alpha(f1,f2,alpha):
    return (1-alpha)*f1+alpha*f2

def lorentzAssymetry(x, mu, sig, alpha):
    model = 1+arctan(alpha*(x-mu)/sig)/(pi/2)
    return model
    
    
def pseudoVoigt(x, mu, sig, A, alpha):
    model1 = alpha * 1/(1 + ((x-mu)/sig)**2) 
    model2 = (1-alpha)*exp(-log(2)*((x-mu)/sig)**2)
    return A* (model1 + model2)

def asymGauss(x, ampl, center, sigma, alpha):
    model1 = exp(-((x-center)**2)/(2*sigma**2))
    model2 = (1+erf((alpha*(x-center))/(sigma*sqrt(2))))
    return ampl * model1 * model2

def benOckoAsymGaussLog(x, qc,a ,b, c, d, sig, ampl, center, sigmaL, sigmaR, scale, offset):
    x = scale*x+offset
    model1 = benOcko(x, qc,a ,b, c, d, sig)
    model2 = asymGauss(x, ampl, center, sigmaL, sigmaR)
    return log10(model1 + model2)

def benOckoAsymSincLog(x, qc,a ,b, c, d, sig, ampl, center, sigmaL, sigmaR, scale, offset):
    x = scale*x+offset
    model1 = benOcko(x, qc,a ,b, c, d, sig)
    model2 = asymSincSqd(x, ampl, center, sigmaL, sigmaR)
    return log10(model1 + model2)

def DyRefAsymGaussLog(x, thickness, roughness, roughnessSub, I0, bgd, resol, energy, relaxation, ampl, center, sigma, alpha, scale, offset):
    x = scale*x+offset
    model1 = DyRefFunction(x, thickness, roughness, roughnessSub, I0, bgd, resol, energy, relaxation)
    model2 = asymGauss(x, ampl, center, sigmaL, sigmaR)
    return log10(model1 + model2)
    
    
    