import sys
sys.path.append('../lib')
import time
import argparse
from recurrence import *
from decomp import *
from quad import *
from sti import *
from itertools import product
from poly import *
from ppoly import *
from fea1d import *
from egrhs import *
from aposteriori import *

def average(x):
    assert len(x) > 0
    return fsum(x) / mpf(len(x))

def pearson(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / sqrt(xdiff2 * ydiff2)


parser = argparse.ArgumentParser(description='Implicit aposteriori benchmark.',add_help=False)
parser.add_argument('-kappa2', metavar='KAPPA2', type=float, required=True, help='Constant in DE.')
parser.add_argument('-nel', metavar='NEL', type=int, required=True, help='Finite elements.')
parser.add_argument('-p', metavar='P', type=int, required=True, help='Finite element degree')
parser.add_argument('-case', metavar='PARAM', required=True, help='Test case name.')
parser.add_argument('-iapnel', metavar='IAPNEL', type=int, required=True, help='IAP elements.')
parser.add_argument('-iapp', metavar='IAPP', type=int, required=True, help='IAP element degree')
args = parser.parse_args()

mp.dps = 20
kappa2 = args.kappa2
a=-1
b=1
nel=args.nel
p=args.p
iapnel=args.iapnel
iapp=args.iapp

bt = 0

if args.case == "spike":
    s = 100
    exact = spike(s)
    dexact = spikeder(s)
    rhs = spikerhs(s,kappa2)
    ad = 0
    bd = 0
    titl = r"$(1-x^2)^{{{0:s}}}$".format(nstr(s))
elif args.case == "arctan":
    r = 100
    exact = arctanjump(r)
    dexact = arctanjumpder(r)
    rhs = arctanjumprhs(r,kappa2)
    ad = exact(a)
    bd = exact(b)
    titl = r'$\frac{{2}}{{\pi}}\arctan({0:s}x)$'.format(nstr(r))
elif args.case == "arctanbub":
    r = 100
    exact = arctanbubjump(r)
    rhs = arctanbubjumprhs(r,kappa2)
    ad = exact(a)
    bd = exact(b)
    titl = r'$\frac{{2}}{{\pi}}\arctan({0:s}x)(1-x^2)$'.format(nstr(r))
elif args.case == "expspike":
    s = 0.01
    exact = expspike(s)
    dexact = lambda x : (-2)/s**2 * x*exp(-x**2/s**2)
    rhs = expspikerhs(s,kappa2)
    ad = exact(a)
    bd = exact(b)
    titl = r"$e^{{-x^2/{0:s}^2}}$".format(nstr(s))
elif args.case == "loggap":
    s = 0.0001
    exact = loggap(s)
    dexact = lambda x : 1.0/(1+x+s)
    rhs = loggaprhs(s,kappa2)
    ad = exact(a)
    bd = exact(b)
    titl = r"$\log(1+x+{0:s})$".format(nstr(s))
elif args.case == "homo":
    ad = 0
    bd = 1
    kappa = sqrt(kappa2)
    c2 = exp(kappa)/(1.0-exp(4.0*kappa))
    c1 = -c2*exp(2*kappa)
    exact = lambda x : c1*exp(kappa*x)+c2*exp(-kappa*x)
    dexact = lambda x: c1*exp(kappa*x)-c2*exp(-kappa*x)
    rhs = lambda phi,r0,r1,el0,el1 : 0.0

    titl = r'$0$'
elif args.case == "sinc":
    d = 5
    exact = lambda x : sin(pi*d*x)/(pi*d*x) if x != 0 else 1.0
    dexact = lambda x : (pi*d*x*cos(pi*d*x)-sin(pi*d*x))/(pi*d*x**2) if x != 0 else 0
    rhs = rhsmpmathquad(lambda x : (2*pi*d*x*cos(pi*d*x)-(2-pi**2*d**2*x**2)*sin(pi*d*x))/(pi*d*x**3) + kappa2*sin(pi*d*x)/(pi*d*x))
    ad = exact(a)
    bd = exact(b)
    titl = "sinc"
else:
    print "error: invalid test case specified. specify one of the following."
    print "spike arctan arctanbub expspike loggap homo"
    quit()

X = linspace(a,b,nel+1)

ab = r_jacobi(p+1)
xw = gauss(ab)

print ">>>> Finite element analysis for p=", p
phi = lagrangecheby(p)
els,G,x,phi = fea1dh(X, lagrangecheby(p), kappa2, bt, [ad,bd], rhs)
uh = ppolyfea1sol(els,G,x,phi)

print ">>>> Implicit a posteriori error estimation"
duh = ppolyder(uh)
dduh = ppolyder(duh)

duX = map(dexact, X)

print ">>>>>> Neumann problem with exact derivatives"
t1 = time.time()
errhhneu,errhhneunorm = iapneu(els,G,x,phi,kappa2,rhs,duX,iapnel,iapp)
t2 = time.time()

print ">>>>>> Neumann problem with face averaging"
t3 = time.time()
errhhneufa,errhhneufanorm = iapneu(els,G,x,phi,kappa2,rhs,None,iapnel,iapp)
t4 = time.time()

print ">>>>>> Dirichlet problem with bubble polynomials"
errhhbub,errhhbubnorm,tbub = iapbub(els,G,x,phi,kappa2,rhs,duX,iapnel,iapp)

errhnorm = ppolyerrh1norm2intv(uh,exact,dexact)

for e in range(len(els)):
    print "{0:>2d} {1:>15s} {2:>15s} {3:>15s} {4:>15s}".format(e, nstr(errhhneunorm[e]), nstr(errhhneunorm[e]), nstr(errhhneufanorm[e]), nstr(errhhbubnorm[e]))

print "{0:f} {1:f} {2:f}".format(t2-t1, t4-t3, tbub)

print "rho(exact,neudexact) =", pearson(errhhneunorm,errhnorm)
print "rho(exact,neufa) =", pearson(errhhneufanorm,errhnorm)
print "rho(exact,bub) =", pearson(errhhbubnorm,errhnorm)


