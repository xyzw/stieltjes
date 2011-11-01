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
        mom0[k] = quadpq(mon(k),mu0,xwl) # int x^k dmu0
        mom1[k] = quadpq(mon(k),mu1,xwl) # int x^k dmu1
    return mom0,mom1

mp.dps = 20

kappa2 = mpf(10)
n = 5
alpha = 1
gln = 2*alpha

print ">>>> Generating Gauss-Legendre quadrature rules for n =", gln
xwl = gauss(r_jacobi(gln))

print "[-1,1]:", quadpq(bw(alpha), ones(1), xwl)
a = 0.0
b = 1.0
print quadapq(bw(alpha), ones(1), a, b, xwl)


xx = linspace(-1,1,1000)
p = bw(alpha)
#print paff
paff = polyaff(p,-1.0,1.0,a,b)
yy = polyvalv(paff,xx)
pl.plot(xx,yy)

pl.show()
