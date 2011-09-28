import sys
sys.path.append('../lib')
from bisection import *
import time

def tridi101(dpsv, nv, epsv):
    print ">>>>>> tridi(1,0,1)"
    print "{0:>4s} {1:>4s} {2:>10s} {3:>15s} {4:>5s} {5:>10s}".\
          format("n", "dps", "eps",  "errmaxnorm", "z", "sec")

    for dpsi in dpsv:
        for ni in nv:
            for epsi in epsv:
                mp.dps = dpsi
                
                lk = lambda k: 2*cos((ni-k+1)*pi/(ni+1))
                se = matrix([map(lk, linspace(1,ni,ni))])

                t1 = time.time()
                s, z = bisection(zeros(ni,1), ones(ni-1,1), epsi, 0, ni-1)
                t2 = time.time()
                err = norm(s-se, inf)
                
                print "{0:4d} {1:4d} {2:>10s} {3:>15s} {4:5d} {5:10f}".format(ni, dpsi, \
                                                                      nstr(epsi), \
                                                                      nstr(err),
                                                                      z, \
                                                                      (t2-t1))


print ">> Eigenproblem"

print ">>>> bisection"

dpsv = [10, 15, 20]
nv = [5, 10, 20]
epsv = [mpf(1e-10), mpf(1e-25), mpf(1e-20)]

tridi101(dpsv, nv, epsv)

