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


g = mpf(100)
n = 5

xw = gauss(r_jacobi(n))

XW = row_join( row_join(xw[:,0],xw[:,0]), row_join(xw[:,1],g*xw[:,1]) )
B = dstis1(n,XW,[n, n],0)

print B

P = r_to_poly(B)

xx = linspace(-1,1,100)
for k in range(n):
    yy = polyvalv(P[k,:],xx)
    pl.plot(xx,yy,label="k="+str(k))

pl.legend()
pl.axis([-1,1,-1.5,1.5])
pl.grid(True)
pl.show()
