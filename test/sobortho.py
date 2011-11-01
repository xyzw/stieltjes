import sys
sys.path.append('../lib')
from recurrence import *
from decomp import *
from quad import *
from sti import *
from chebyshev import *
from itertools import product
from poly import *
from ppoly import *
from fea1d import *
import pylab as pl
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':'18'})

# Bubble weight (1-x^(2*alpha))
def bw(alpha=1):
    w = -mon(2*alpha)
    w[2*alpha] = mpf(1)
    return w

# Calculate modified moments
def bwmom(n,kappa2,alpha,xwl):
    w = bw(alpha)
    mu0 = polyaxpy(kappa2,conv(w,w),conv(polyder(polyder(w)),-w)) # k^2w^2-w''w
    mu1 = conv(w,w)
    mom0 = zeros(1,2*n);
    mom1 = zeros(1,2*n);
    for k in range(0,2*n):
        mom0[k] = quadpq(mon(k),mu0,xwl)
        mom1[k] = quadpq(mon(k),mu1,xwl)
    return mom0,mom1

mp.dps = 20

kappa2 = mpf(10)
n = 5
alpha = 1
gln = n+2*alpha

mp.dps = 50
print ">>>> Generating Gauss-Legendre quadrature rules for n =", gln
xwl = gauss(r_jacobi(gln))
mp.dps = 20

print ">>>> Generating Bubble polynomials for kappa2={0:s}, n={1:d}, alpha={2:d}".format(nstr(kappa2),n,alpha)
mom0,mom1 = bwmom(n, kappa2, alpha, xwl)
ab0,n2 = chebyshev(n,mom0)
ab1,n2 = chebyshev(n,mom1)
xw0 = gauss(ab0)
xw1 = gauss(ab1)

XW = row_join( row_join(xw0[:,0],xw1[:,0]), row_join(xw0[:,1],xw1[:,1]) )
B,n2 = dstis1(n,XW,[n, n],0)
P = r_to_poly(B)

P = convm(P, bw(alpha))

# H1-normalize
for k in range(n):
    h1ni = mpf(1)/sqrt(h1norm2(P[k,:],xwl))
    P[k,:] *= h1ni

G = zeros(n)
for k in range(n):
    for l in range(n):
        G[k,l] = h1inner(P[k,:],P[l,:],xwl,kappa2)
     
print chop(G)

xx = linspace(-1,1,1000)
for k in range(n):
    yy = polyvalv(P[k,:],xx)
    pl.plot(xx,yy,label=r"$k={0:d}$".format(k),linewidth=2.0)

pl.title(r"$\{{C_k\}}$ bubor\'ekpolinomok ($\kappa^2={0:s}$, $n={1:d}$, $\alpha={2:d})$".format(nstr(kappa2),n,alpha))
pl.legend()
pl.axis([-1,1,-0.6,0.6])
pl.grid(True)
pl.show()
