import sys
sys.path.append('../lib')
from recurrence import *
from decomp import *
from quad import *
from sti import *
from itertools import product
from poly import *
from fea1d import *

import time

def testjac(dpsv, nv, av, bv, epsv):
    print ">>>>>> r_jacobi(n,a,b)"
    print "{0:>4s} {1:>4s} {2:>4s} {3:>4s} {4:>10s} {5:>15s} {6:>10s} {7:>10s}".\
          format("n", "a", "b", "dps", "eps",  "errmaxnorm", "tg", "ts")

    for ni,ai,bi,dpsi,epsi in product(nv,av,bv,dpsv,epsv):
        mp.dps = dpsi
        
        # Compute using the closed form
        ab = r_jacobi(ni,ai,bi)
        
        # Calculate Gauss quadrature rules
        t1 = time.time()
        xw = gauss(ab)
        t2 = time.time()
        
        # Invoke the Stieltjes algorithm to calculate ab again
        t3 = time.time()
        ab1 = dsti(xw)
        t4 = time.time()
        
        err = norm(ab-ab1, inf)
        
        print "{0:>4d} {1:>4s} {2:>4s} {3:>4d} {4:>10s} {5:>15s}  {6:>10f} {7:>10f}".format(\
            ni, nstr(ai), nstr(bi), dpsi, nstr(epsi), nstr(err), t2-t1, t4-t3)

print ">> Gauss quadrature rule generation"

print ">>>> gauss"

dpsv = [5, 15, 50]
nv = [10, 20]
av = [mpf(1)]
bv = [mpf(1)]
epsv = [mpf(1e-5), mpf(1e-15), mpf(1e-50)]

#testjac(dpsv,nv,av,bv,epsv)

n=5
kappa = 1
ab = r_jacobi(n,0,0)
p = ab_to_poly(ab)
xw = gauss(ab)

els = uniform_els(-1,1,3)
print els
fea_diri2(els, kappa, 0, xw)

quit()

from pylab import *
xx = arange(-1.0, 1.0, 0.01)
for i in range(n):
    plot(xx,polyvalv(p[i,:],xx))
axis([-1,1,-0.5,0.5])
grid(True)
show()

#

#XW = row_join( row_join(xw[:,0],xw[:,0]), row_join(xw[:,1],xw[:,1]) )
#B = pdstis1(XW,[n, n],1)
#print B


