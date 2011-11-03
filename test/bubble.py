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
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':'18'})

                
def testprebub(nv, alphav, hv, kappa2v):
    print ">>>>>> prebub"
    print "{0:>2s} {1:>2s} {2:>8s} {3:>8s} {4:>10s} {5:>6s} {6:>6s}".\
          format("n", "alpha", "h", "kappa2", "errmax", "tgl", "tpb")

    mp.dps = 50
    xwll = gauss(r_jacobi(max(nv)+2*max(alphav)))
    mp.dps = 20
    errmax = mpf(-inf)

    for ni,alphai,hi,kappa2i, in product(nv,alphav,hv,kappa2v):

        mp.dps = 50
        gln = ni+2*alphai
        
        t1 = time.time()
        xwl = gauss(r_jacobi(gln))
        t2 = time.time()

        mp.dps = 20

        t3 = time.time()
        P = prebub(ni,hi,kappa2i,alphai,xwl)
        t4 = time.time()

        a = -mpf(0.5)*hi
        b = mpf(0.5)*hi
        G = gramh1(P,a,b,kappa2i,xwll)
        err = maxod(G)
        errmax = max(errmax,err)
        print "{0:>2d} {1:>5d} {2:>8s} {3:>8s} {4:>10s} {5:>6.3f} {5:>6.3f}".format(\
            ni, alphai, nstr(hi), nstr(kappa2i), nstr(err,4), (t2-t1)*1000.0, (t4-t3)*1000.0)
        #print chop(G)

    print ">> errmax = {0:s}".format(nstr(errmax))

nv = [2, 5]
alphav = [1]
hv = [mpf(2), mpf(10), mpf(0.01), mpf(0.0001), mpf(1000)]
kappa2v = [mpf(1), mpf(10), mpf(100), mpf(0.0001)]
testprebub(nv,alphav,hv,kappa2v)

## xx = linspace(a,b,1000)
## for k in range(n):
##     yy = polyvalv(P[k,:],xx)
##     pl.plot(xx,yy,label=r"$k={0:d}$".format(k),linewidth=2.0)

## pl.title(r"$\{{C_k\}}$ bubor\'ekpolinomok ($\kappa^2={0:s}$, $n={1:d}$, $\alpha={2:d})$".format(nstr(kappa2),n,alpha))
## pl.legend()
## pl.axis([float(a),float(b),-1,1])
## pl.grid(True)
## pl.show()
