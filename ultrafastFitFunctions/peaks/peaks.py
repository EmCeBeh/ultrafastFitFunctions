from numpy import *  # pi, exp, sqrt, etc.
from scipy.special import erf, erfc

from ..backgrounds import benOcko1, benOcko2


def lorentz(x, mu, sig, A, c):
    model = A / (1 + ((x - mu) / sig) ** 2)
    return model


def alpha(f1, f2, alpha):
    return (1 - alpha) * f1 + alpha * f2


def lorentzAssymetry(x, mu, sig, alpha):
    model = 1 + arctan(alpha * (x - mu) / sig) / (pi / 2)
    return model


def pseudoVoigt(x, mu, sig, A, alpha):
    model1 = alpha * 1 / (1 + ((x - mu) / sig) ** 2)
    model2 = (1 - alpha) * exp(-log(2) * ((x - mu) / sig) ** 2)
    return A * (model1 + model2)


def asymGauss(x, ampl, center, sigma, alpha):
    model1 = exp(-((x - center) ** 2) / (2 * sigma**2))
    model2 = 1 + erf((alpha * (x - center)) / (sigma * sqrt(2)))
    return ampl * model1 * model2


def benOckoAsymGaussLog(x, qc, a, b, c, d, sig, ampl, center, sigmaL, sigmaR):
    model1 = benOcko(x, qc, a, b, c, d, sig)
    model2 = asymGauss(x, ampl, center, sigmaL, sigmaR)
    return log10(model1 + model2)


def benOckoAsymSincLog1(x, qc, a, b, c, d, sig, ampl, center, sigmaL, sigmaR):
    """
    qc is where the background starts
    """
    model1 = benOcko1(x, qc, a, b, c, d, sig)
    model2 = asymSincSqrd(x, ampl, center, sigmaL, sigmaR)
    return log10(model1 + model2)


def benOckoAsymSincLog2(x, qc, a, b, c, d, sig, ampl, center, sigmaL, sigmaR):
    """
    qc is where the background starts
    """
    model1 = benOcko2(x, qc, a, b, c, d, sig)
    model2 = asymSincSqrd(x, ampl, center, sigmaL, sigmaR)
    return log10(model1 + model2)


def DyRefAsymGaussLog(
    x,
    thickness,
    roughness,
    roughnessSub,
    I0,
    bgd,
    resol,
    energy,
    relaxation,
    ampl,
    center,
    sigma,
    alpha,
    scale,
    offset,
):
    x = scale * x + offset
    model1 = DyRefFunction(
        x, thickness, roughness, roughnessSub, I0, bgd, resol, energy, relaxation
    )
    model2 = asymGauss(x, ampl, center, sigmaL, sigmaR)
    return log10(model1 + model2)


def asymSincSqrd(x, ampl, center, sigmaL, sigmaR):
    sincL = sinc((x[x < center] - center) / sigmaL) ** 2
    sincR = sinc((x[x >= center] - center) / sigmaR) ** 2
    return ampl * r_[sincL, sincR]


def sincSqrd(x, ampl, center, sigma):
    return ampl * sinc((x - center) / sigma) ** 2
