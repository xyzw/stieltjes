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

                
def testprebub(nv, alphav, hv, kappa2v, pv, nelv):
    print ">>>>>> prebub"
    print "{0:>4s} {1:>10s} {2:>4s} {3:>10s} {4:>8s} {5:>8s} {6:>10s}".\
          format("n", "alpha", "p", "nel", "h", "kappa2", "errmax")

    mp.dps = 50
    xwll = gauss(r_jacobi(max(nv)+2*max(alphav)))
    mp.dps = 20
    errmax = mpf(-inf)

    for ni,alphai,hi,kappa2i,pi,neli in product(nv,alphav,hv,kappa2v,pv,nelv):

        mp.dps = 50
        gln = ni+2*alphai
        
        xwl = gauss(r_jacobi(gln))

        mp.dps = 20

        D0,D0inner = bubblexnortho(ni,hi,kappa2i,alphai,pi,neli)
        D1,D1inner = bubbleortho(ni,hi,kappa2i,alphai,pi,neli)

        G0 = ppolygramh1(D0,kappa2i,xwll)
        G1 = ppolygramh1(D1,kappa2i,xwll)

        print ">>>> G0"
        print chop(G0)
        print ">>>> G1"
        print chop(G1)

        err0 = maxod(G0)
        err1 = maxod(G1)
        
        errmax = max(errmax,err0)
        print "{0:>2d} & {1:>5d} & {2:>5d} & {3:>5d} & {4:>8s} & {5:>8s} & \\verb!{6:>10s}! & \\verb!{7:>10s}! \\\\".format(\
            ni, alphai, pi, neli, nstr(hi), nstr(kappa2i), nstr(err0,4), nstr(err1,4))
        #print chop(G)

    print ">> errmax = {0:s}".format(nstr(errmax))

nv = [5, 10]
alphav = [10, 5]
hv = [ mpf(0.0001), mpf(1000)]
kappa2v = [ mpf(1)]
pv = [1, 2]
nelv = [2, 3]
testprebub(nv,alphav,hv,kappa2v,pv,nelv)
