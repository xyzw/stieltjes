import sys
sys.path.append('../lib')
from recurrence import *
from decomp import *
from gauss import *
from sti import *
from itertools import product
from poly import *
from fea1d import *
import pylab as pl
from matplotlib.widgets import Slider


n=5
kappa = 1
ab = r_jacobi(n,0,0)
p = ab_to_poly(ab)
xw = gauss(ab)

els = equidistant(-1,1,3)
print els
fea_diri2(els, 2, kappa, matrix([-1,0,1]), xw)

print ab_to_poly(ab)

## f = matrix([-1,0,1])
## ftilde = polyaff(f, -0.5, 0.5, -1, 1)
## xx = pl.arange(-1.0, 1.0, 0.01)
## pl.plot(xx,polyvalv(f,xx))
## pl.plot(xx,polyvalv(ftilde,xx))
## pl.axis([-1,1,-1,1])
## pl.grid(True)
## pl.show()

#quit()

ax = pl.subplot(111)
pl.subplots_adjust(bottom=0.25)

pl.axis([-1,1,-2,2])
pl.grid(True)
ps = Slider(pl.axes([0.1, 0.1, 0.7, 0.03]), 'n', 2, 15, valinit=2, valfmt='%1.0f')

def update(val):
    ax.clear()
    p = int(ps.val)
    L = lobatto(p)
    xx = pl.linspace(-1, 1, 100)
    ax.plot(pl.linspace(-1,1,p),pl.zeros(p),'o')
    ax.plot(pl.linspace(-1,1,p),pl.ones(p),'o')
    for i in range(p):
        ax.plot(xx,polyvalv(L[i,:],xx))
    pl.draw()
    
ps.on_changed(update)

pl.title('Lobatto basis polynomials')
pl.show()

