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

mp.dps = 15

kappa2 = 1
a=-1
b=1
bt = 0
s = 100
exact = arctanjump(s)
dexact = arctanjumpder(s)
rhs = arctanjumprhs(s,kappa2)
ad = exact(a)
bd = exact(b)

#d = 10
#exact = lambda x : sin(pi*d*x)/(pi*d*x) if x != 0 else 1.0
#dexact = lambda x : (pi*d*x*cos(pi*d*x)-sin(pi*d*x))/(pi*d*x**2) if x != 0 else 0
#rhs = rhsmpmathquad(lambda x : (2*pi*d*x*cos(pi*d*x)-(2-pi**2*d**2*x**2)*sin(pi*d*x))/(pi*d*x**3) + kappa2*sin(pi*d*x)/(pi*d*x) if x != 0 else kappa2)
#ad = exact(a)
#bd = exact(b)

print "{0:>6s} & {1:>6s} & {2:>6s} & {3:>6s} & {4:>25s} & {5:>25s} & {6:>25s} & {7:>25s} & {8:>10s} & {9:>10s} \\\\".format(\
    "nel","p","iapnel","iapp","rhoneufa","rhobub","rhobuba","rhoxnbub","dtneufa","dtbub")

for neli in range(10,20):
    for pi in range(1,2):
        for iapneli in range(3,4):
            for iappi in range(3,4):

                X = linspace(a,b,neli+1)

                ab = r_jacobi(pi+1)
                xw = gauss(ab)

                #print ">>>> Finite element analysis for p=", p
                phi = lagrangecheby(pi)
                els,G,x,phi = fea1dh(X, lagrangecheby(pi), kappa2, bt, [ad,bd], rhs)
                uh = ppolyfea1sol(els,G,x,phi)

                #print ">>>> Implicit a posteriori error estimation"
                duh = ppolyder(uh)
                dduh = ppolyder(duh)

                #print ">>>>>> Neumann problem with exact derivatives"
                #t1 = time.time()
                #errhhneu,errhhneunorm = iapneu(els,G,x,phi,kappa2,rhs,duX,iapneli,iappi)
                #t2 = time.time()

                #print ">>>>>> Neumann problem with face averaging"
                t3 = time.time()
                errhhneufa,errhhneufanorm = iapneu(els,G,x,phi,kappa2,rhs,None,iapneli,iappi)
                t4 = time.time()

                #print ">>>>>> Dirichlet problem with bubble polynomials"
                #errhhbub,errhhbubnorm,tbub = iapbub(els,G,x,phi,kappa2,rhs,iapneli,iappi)

#                errhhbubxn,errhhbubxnnorm,tbubxn = iapbub(els,G,x,phi,kappa2,rhs,iapneli,iappi,False)

#                errhhbubxnspread,errhhbubxnspreadnorm,tbubxnspread = iapbubspread(els,G,x,phi,kappa2,rhs,iapneli,iappi,False)

                errhhbubortho,errhhbuborthonorm,tbubortho = iapbubortho(els,G,x,phi,kappa2,rhs,iapneli,iappi)
                errhhbuba,errhhbubanorm,tbuba = iapbuba(els,G,x,phi,kappa2,rhs,iappi)
                errhhbubxnortho,errhhbubxnorthonorm,tbubxnortho = iapbubxnortho(els,G,x,phi,kappa2,rhs,iapneli,iappi)
                                
                errhnorm = ppolyerrh1norm2intv(uh,exact,dexact)

                #for e in range(len(els)):
                #    print "{0:>2d} {1:>15s} {2:>15s} {3:>15s} {4:>15s}".format(e, nstr(errhhneunorm[e]), nstr(errhhneunorm[e]), nstr(errhhneufanorm[e]), nstr(errhhbubnorm[e]))

                print "{0:>6d} & {1:>6d} & {2:>6d} & {3:>6d} & {4:>25s} & {5:>25s} & {6:>25s} & {7:>25s} & {8:>10f} & {9:>10f} \\\\".\
                      format(neli,pi,iapneli,iappi, pearson(errhhneufanorm,errhnorm), pearson(errhhbuborthonorm,errhnorm),pearson(errhhbubanorm,errhnorm), pearson(errhhbubxnorthonorm,errhnorm),\
                             t4-t3, tbubortho)

                ## print "{0:>6d} & {1:>6d} & {2:>6d} & {3:>6d} & {4:>25s} & {5:>25s} & {6:>25s} & {7:>25s} & {8:>25s} & {9:>10f} & {10:>10f} & {11:>10f} & {12:>10f} & {13:>10f}\\\\".\
                ##       format(neli,pi,iapneli,iappi,\
                ## pearson(errhhneufanorm,errhnorm),\
                ## pearson(errhhbubnorm,errhnorm), \
                ## pearson(errhhbubxnnorm,errhnorm), \
                ## pearson(errhhbubxnspreadnorm,errhnorm),\
                ## pearson(errhhbuborthonorm,errhnorm), t4-t3, tbub, tbubxn, tbubxnspread, tbubortho)

                #print "{0:f} {1:f} {2:f}".format(t2-t1, t4-t3, tbub)

                #print "rho(exact,neudexact) =", pearson(errhhneunorm,errhnorm)
                #print "rho(exact,neufa) =", pearson(errhhneufanorm,errhnorm)
                #print "rho(exact,bub) =", pearson(errhhbubnorm,errhnorm)


