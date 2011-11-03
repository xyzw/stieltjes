import sys
sys.path.append('../lib')
import time
from recurrence import *
from decomp import *
from quad import *
from sti import *
from chebyshev import *
from itertools import product
from poly import *
from ppoly import *
from bub import *
from fea1d import *
import pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
#rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':'18'})

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

def rhsbubble(C, h, kappa2, xw):
    f = polyaxpy(mpf(-1), polyder(polyder(C)), kappa2*C) # -C'' + kappa2*C
    return lambda phi,el0,el1 : ((el1-el0)/mpf(2)) * quadpq(polyaff(f, mpf(-1), mpf(1), el0, el1),phi,xw)


pl.grid(True)

mp.dps = 20

# Generate Bubble polynomials
kappa2=mpf(100)
h=mpf(10)
n=5
nel=3
alpha = 1
r=500
p=3

xwl = gauss(r_jacobi(n+2*alpha))
P = prebub(n,h,kappa2,alpha,xwl)

xwr = gauss(r_jacobi(n+2*alpha+2))

X = linspace(-mpf(0.5)*h,mpf(0.5)*h,nel+1)

for k in range(n):
    print k

    els,G,x,phi = fea_diri2(X, p, kappa2, [mpf(0),mpf(0)], rhsbubble(P[k,:], h, kappa2, xwl), xwr)
    pp = ppolyfea1sol(els,G,x,phi)
    #pf = polytoppoly(poly([-1,0,1]),els,a,b)
    #perr = ppolyaxpy(-1,pp,pf)
    #ploterr(p,[perr])
    #print "l2err =", ppolyl2norm(perr,xw)

    Chat = ppolyaxpy(-1,polytoppoly(P[k,:],pp.intv,-mpf(0.5)*h,mpf(0.5)*h),pp)

    xx,yy = ppolyval(Chat,50)
    
    #pl.plot(xx,polyvalv(P[k,:],xx),label=r"$C_{0:d}$".format(k),linewidth=0.8)
    pl.plot(xx,yy,label=r"$\widehat{{C}}_{0:d}$".format(k),linewidth=1.2)



#rhs = arctanjumprhs(r,kappa2)
#ad = arctanjump(r)(a)
#bd = arctanjump(r)(b)
    
#els,G,x,phi = fea_diri2(X, p, kappa2, [ad,bd], rhs, xw)

#s = 10

#xx = linspace(a,b,1000)

#yy = map(polyspike(s), xx)
#pl.plot(xx,yy,linewidth=2,label="exact")

#yy = map(arctanjump(r), xx)
#pl.plot(xx,yy,linewidth=2,label="exact")


#pl.title(r"Polynomial spike function $(1-x^2)^s$ ($n=" + str(n) + "$, $s=" + str(s) + "$)")

pl.title(r"$\{{\widehat{{C}}_k\}}$ ($\kappa^2={0:s}$, $n={1:d}$, $\alpha={2:d}$, $p={3:d}$)".format(nstr(kappa2),n,alpha,p))
pl.legend()
plotels(X)
#print phi


pl.axis([float(-0.5*h),float(0.5*h),-0.2,0.2])

pl.show()
quit()
