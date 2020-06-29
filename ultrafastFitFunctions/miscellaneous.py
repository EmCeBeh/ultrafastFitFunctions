from numpy import *

class Misc():
    
    def saturation(x, mu, A=1, k=1):
        #if mu < min(x) or mu > max(x):
        #   raise ValueError()
        
        s = sigmoid(x, mu, k=k)
        slope = max((k*s*(1-s)))
        y =  (x-mu)*slope+1/2
        s[x<mu] = y[x<mu]
        
        c = 1/2-mu*slope
        s = s-c
        
        s[x<0] = 0
        s = A*s/max(s)
        
        return s