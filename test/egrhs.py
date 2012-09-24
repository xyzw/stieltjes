from fea1d import *
from poly import *

def spike(s):
    spike = lambda x : (1-x**2)**s
    return spike

def spikeder(s):
    spikeder = lambda x : s*(1-x**2)**(s-1)*(-2*x)
    return spikeder

def spikerhs(s,kappa2):
    return rhsmpmathquad(lambda x : 2*s*(1-x**2)**(s-1)-s*(s-1)*4*x**2*(1-x**2)**(s-2) + kappa2*(1-x**2)**s)

def neuspike(s):
    spike = lambda x : (1-x**2)**s+1
    return spike

def neuspikeder(s):
    spikeder = lambda x : s*(1-x**2)**(s-1)*(-2*x)
    return spikeder

def neuspikerhs(s,kappa2):
    return rhsmpmathquad(lambda x : -2*s*(1-x**2)**(s-2)*((2*s-1)*x**2-1) + kappa2*((1-x**2)**s + 1))

def arctanjump(r):
    atanjump = lambda x : 2/pi*atan(x*mpf(r))
    return atanjump

def arctanjumpder(r):
    atanjump = lambda x : (2/pi)*r/(r**2*x**2+1)
    return atanjump

def arctanjumprhs(r,kappa2):
    atanjumprhs = lambda x : mpf(4)*r**3*x/(pi*(r**2*x**2+1)**2) + kappa2*2/pi*atan(x*r)
    return rhsmpmathquad(atanjumprhs)

def arctanbubjump(r):
    atanjump = lambda x : 2/pi*atan(x*mpf(r))*(1-x**2)
    return atanjump

def arctanbubjumpder(r):
    atanjump = lambda x : 2/pi*((r-r*x**2)/(r**2*x**2+1) -2*x*atan(x*r))
    return atanjump

def arctanbubjumprhs(r,kappa2):
    atanjumprhs = lambda x : (mpf(4)*r*x*(r**2*(x**2+1)+2) + (r**2*x**2+1)**2*atan(r*x))/(pi*(r**2*x**2+1)**2) + kappa2*mpf(2)/pi*atan(x*r)*(mpf(1)-x**2)
    return rhsmpmathquad(atanjumprhs)

def expspike(a):
    return lambda x : exp(-x**2/a**2)

def expspikeder(a):
    return lambda x : exp(-x**2/a**2)*(-2*x/a**2)

def expspikerhs(a,kappa2):
    return rhsmpmathquad(lambda x : (2.0/a**4)*exp(-x**2/a**2)*(a**2-2*x**2) + kappa2*exp(-x**2/a**2))

def loggap(s):
    return lambda x : log(1+x+s)

def loggapder(s):
    return lambda x : 1.0/(1.0+x+s)

def loggaprhs(s,kappa2):
    return rhsmpmathquad(lambda x : 1.0/(1+x+s)**2 + kappa2*log(1+x+s))
