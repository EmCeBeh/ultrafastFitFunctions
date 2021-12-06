from numpy import exp, sqrt, ones_like, pi
from scipy.special import erf, erfc


def step(x):
    #nice name, but renaming to step would make many things shorter
    output = ones_like(x)
    output[x < 0] = 0
    #output[x >= 0] = 1
    return output

def expDecay(x, tau, A):
    return A*step(x)*(exp(-x/tau)-1) + 1
    
def sglExpDecay(x, mu, tau, A, B):
    return A-(1-B*exp(-(x-mu)/tau))*step(mu-x)
    
def expDecayOffs(x, tau, A, c):
    return A*step(x)*(exp(-x/tau)-1) + 1 + c
    
def doubleExpDecay(x, tau1, tau2, A, B):
    model1 = A*step(x)*(exp(-x/tau1)-1)
    model2 = B*step(x)*(exp(-x/tau2)-1)
    return model1 + model2 + 1

def doubleExpDecay_II(x, tau1, tau2, A, B, mu, c):
    return A*step(x-mu)*(exp(-(x-mu)/tau1)-1)+B*step(x-mu)*(exp(-(x-mu)/tau2)-1) + c

def doubleDecay(x, mu, tau1, tau2, A, q):
    B = (1-q)*A    
    A = q*A
    return doubleExpDecay(x-mu,tau1,tau2,A,B)
    
def doubleDecay2(x, mu, tau1, tau2, A, B):
    return doubleExpDecay(x-mu,tau1,tau2,A,B)


def expConvGauss(x, tau, A, sig):
    model1 = exp(-x/tau)*exp(sig**2/(2*tau**2))*(erf((sig**2-x*tau)/(sqrt(2)*sig*tau))-1)
    return -A/2*(model1)

def expConvGaussNormalised(x, tau, sig):
    exp_enu = sig**2-2*(x)*tau
    exp_den = 2*tau**2
    model1  = exp(exp_enu/exp_den)
    erf_enu = sig**2+(-x)*tau
    erf_den = sqrt(2)*sig*tau
    model2  = erf(erf_enu/erf_den)-1
    return 1/2 * model1 * model2

def expConvGauss2(x, mu, tau, A, sig):
    model1 = exp(-(x-mu)/tau)*exp(sig**2/(2*tau**2))*(erf((sig**2-(x-mu)*tau)/(sqrt(2)*sig*tau))-1)
    return -A/2*model1
    
def AHConvGauss(x, A, sig):
    model1 = erf(x/(sqrt(2)*sig))
    return -A/2*(model1-2/A+1)
    
def ABHConvGauss(x, A, B, sig):
    model1 = -A/2*erf(x/(sqrt(2)*sig))
    model2 = -B/2*erf(x/(sqrt(2)*sig))
    return model1+model2+1-A/2-B/2

def expConvGaussApprox(x, tau, A, sig):
    C = sqrt(2/pi)*exp(-x**2/(2*sig**2))
    r = sig/tau
    k = C/r - x*C/(sig*r**2) + C*(sig**2-x**2)/(sig**2*r**3)#+C*x*(3*sig**2-x**2)/(sig**3*r**4)-C*(3*sig**4-6*sig**2*x**2+x**4)/(sig**4*r**5)
    return -A/2*k
    
def doubleDecaySingleConv(x, t0, tau1, tau2, A, B, sig, c):
    C = -(A+B)
    model1 = -A*expConvGaussNormalised(x-t0, tau1, sig)
    model2 = -B*expConvGaussNormalised(x-t0, tau2, sig)
    model3 = +1/2*C*erfc(-(x-t0)/(sqrt(2*sig**2))) +1
    return  model1 + model2 + model3 + c

def doubleDecayDoubleConv(x, mu, tau1, tau2, A, q, alpha, sigS, sigH):

    # A = overall amplitude
    # q = tau1 fraction of A
    # alpha = slicing fraction

    B = (1-q)*A    
    A = q*A

    model1 = expConvGauss(x-mu,tau1,A,sigS)
    model2 = expConvGauss(x-mu,tau2,B,sigS)
    model3 = ABHConvGauss(x-mu,A,B,sigS)
    if sigH/tau1 < 5:
        model4 = expConvGauss(x-mu,tau1,A,sigH)
    else:
        model4 = expConvGaussApprox(x-mu,tau1,A,sigH)

    if sigH/tau2 < 5:
        model5 = expConvGauss(x-mu,tau2,B,sigH)
    else:
        model5 = expConvGaussApprox(x-mu,tau2,B,sigH)
    
    model6 = ABHConvGauss(x-mu,A,B,sigH)

    return alpha*(model1 + model2 + model3)+(1-alpha)*(model4 + model5 + model6)

def doubleDecayDoubleConv2(x, mu, tau1, tau2, A, B, alpha, sigS, sigH):
    
    model1 = expConvGauss(x-mu,tau1,A,sigS)
    model2 = expConvGauss(x-mu,tau2,B,sigS)
    model3 = ABHConvGauss(x-mu,A,B,sigS)
    
    if alpha == 1:
        return model1 + model2 + model3
        
    else:
        
        if sigH/tau1 < 5:
            model4 = expConvGauss(x-mu,tau1,A,sigH)
        else:
            model4 = expConvGaussApprox(x-mu,tau1,A,sigH)

        if sigH/tau2 < 5:
            model5 = expConvGauss(x-mu,tau2,B,sigH)
        else:
            model5 = expConvGaussApprox(x-mu,tau2,B,sigH)
        
        model6 = ABHConvGauss(x-mu,A,B,sigH)

        return alpha*(model1 + model2 + model3)+(1-alpha)*(model4 + model5 + model6)

def doubleDecayConvScale(x, mu, tau1, tau2, A, q, alpha, sigS, sigH, I0):
    B = (1-q)*A    
    A = q*A

    model1 = expConvGauss(x-mu,tau1,A,sigS)
    model2 = expConvGauss(x-mu,tau2,B,sigS)
    model3 = ABHConvGauss(x-mu,A,B,sigS)
    if sigH/tau1 < 5:
        model4 = expConvGauss(x-mu,tau1,A,sigH)
    else:
        model4 = expConvGaussApprox(x-mu,tau1,A,sigH)

    if sigH/tau2 < 5:
        model5 = expConvGauss(x-mu,tau2,B,sigH)
    else:
        model5 = expConvGaussApprox(x-mu,tau2,B,sigH)
    
    model6 = ABHConvGauss(x-mu,A,B,sigH)

    return I0*(alpha*(model1 + model2 + model3)+(1-alpha)*(model4 + model5 + model6))

def DecayConv(x, mu, tau, A, sig):
    
    model1 = expConvGauss(x-mu,tau,A,sig)
  
    model3 = AHConvGauss(x-mu,A,sig)
    
    return model1 + model3

def doubleDecayConvSqrd(x, mu, tau1, tau2, A, q, alpha, sigS, sigH):
    B = (1-q)*A    
    A = q*A

    model1 = expConvGauss(x-mu,tau1,A,sigS)
    model2 = expConvGauss(x-mu,tau2,B,sigS)
    model3 = ABHConvGauss(x-mu,A,B,sigS)
    if sigH/tau1 < 5:
        model4 = expConvGauss(x-mu,tau1,A,sigH)
    else:
        model4 = expConvGaussApprox(x-mu,tau1,A,sigH)

    if sigH/tau2 < 5:
        model5 = expConvGauss(x-mu,tau2,B,sigH)
    else:
        model5 = expConvGaussApprox(x-mu,tau2,B,sigH)
    
    model6 = ABHConvGauss(x-mu,A,B,sigH)

    return (alpha*(model1 + model2 + model3)+(1-alpha)*(model4 + model5 + model6))**2
