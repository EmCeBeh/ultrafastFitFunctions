from numpy import exp

def sigmoid(x, mu, A, k, c):
    model = A/(1+exp(-k*(x-mu)))+c
    return model

def saturation(x, mu, A=1, k=1):
    #if mu < min(x) or mu > max(x):
    #   raise ValueError()
    
    s = sigmoid(x, mu, 1, k, 0)
    slope = max((k*s*(1-s)))
    y =  (x-mu)*slope+1/2
    s[x<mu] = y[x<mu]
    
    c = 1/2-mu*slope
    s = s-c
    
    s[x<0] = 0
    s = A*s/max(s)
    return s

def saturation2(x, mu, A=1, k=1):
    
    #if mu < min(x) or mu > max(x):
    #   raise ValueError()
    
    s = sigmoid(x, mu, A=1,k=k, c=0)
    slope = k/4
    
    y =  (x-mu)*slope+1/2
    s[x<mu] = y[x<mu]
    
    c = 1/2-mu*slope
    s = s-c
    s[x<0] = 0
    s = A*s
    
    return s