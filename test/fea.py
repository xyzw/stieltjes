import sys
sys.path.append('../lib')
from recurrence import *
from decomp import *
from quad import *
from sti import *
from itertools import product
from poly import *
from ppoly import *
from fea1d import *
import pylab as pl
from matplotlib.widgets import Slider
from matplotlib.pyplot import legend
from egrhs import *


# Plot element bars
def plotels(X):
    for l in range(len(X)):
        pl.plot([X[l],X[l]],[-1000,1000],color='grey',linestyle='--')
        
# Plot errors
def ploterr(p,pp):
    for q in pp:
        xx,yy = ppolyval(q,50)
        pl.plot(xx,yy,label="p="+str(p))

# Plot basis elements
def plotbasis(phi):
    xx = linspace(-1,1,100)
    for i in range(len(phi)):
        pl.plot(xx,polyvalv(phi[i,:],xx))

pl.grid(True)

mp.dps = 40

kappa2 = 1
a=-1
b=1
n=3
r=100
s=10

X = linspace(a,b,n+1)

#rhs = arctanjumprhs(r,kappa2)
#ad = arctanjump(r)(a)
#bd = arctanjump(r)(b)

testcase = sys.argv[1]

if testcase == "spike":
    s = 10
    exact = spike(s)
    rhs = spikerhs(s,kappa2)
    ad = 0
    bd = 0
elif testcase == "arctan":
    r = 100
    exact = arctanjump(r)
    rhs = arctanjumprhs(r,kappa2)
    ad = exact(a)
    bd = exact(b)
elif testcase == "arctanbub":
    r = 100
    exact = arctanbubjump(r)
    rhs = arctanbubjumprhs(r,kappa2)
    ad = exact(a)
    bd = exact(b)
elif testcase == "expspike":
    s = 0.01
    exact = expspike(s)
    rhs = expspikerhs(s,kappa2)
    ad = exact(a)
    bd = exact(b)
else:
    print "error: invalid test case specified. specify one of the following."
    print "spike arctan arctanbub"
    quit()

#rhs = rhsmpmathquad(lambda x : 3-x**2)


for p in range(1,10):
    ab = r_jacobi(p+1)
    xw = gauss(ab)
    print p

    phi, phid = lagrangecheby(p)
    
    els,G,x,phi = fea_diri2(X, phi, phid, kappa2, 0, [ad,bd], rhs, xw)

    pp = ppolyfea1sol(els,G,x,phi)
 #   pf = polytoppoly(poly([-1,0,1]),els,a,b)
#    perr = ppolyaxpy(-1,pp,pf)
    #ploterr(p,[perr])
    #print "l2err =", ppolyl2norm(perr,xw)

    xx,yy = ppolyval(pp,50)
    pl.plot(xx,yy,label="p=" + str(p))


xx = linspace(a,b,1000)
yy = map(exact, xx)
pl.plot(xx,yy,linewidth=2,label="exact")

pl.legend()

plotels(X)
#print phi
#plotbasis(phi)

pl.axis([-1,1,-2,2])

pl.show()

