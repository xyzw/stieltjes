import sys
sys.path.append('../lib')
from recurrence import *
from decomp import *
from gauss import *
from sti import *
from itertools import product
from poly import *
from ppoly import *
from fea1d import *
import pylab as pl
from matplotlib.widgets import Slider

# Plot element bars
def plotels(X):
    for l in range(len(X)):
        pl.plot([X[l],X[l]],[-1000,1000],color='grey',linestyle='--')
        
# Plot errors
def ploterr(pp):
    for p in pp:
        xx,yy = ppolyval(p,50)
        pl.plot(xx,yy)

# Plot basis elements
def plotbasis(phi):
    xx = linspace(-1,1,100)
    for i in range(len(phi)):
        pl.plot(xx,polyvalv(phi[i,:],xx))


pl.grid(True)

mp.dps = 20

for p in range(1,8):
    a=-1
    b=1
    n=2
    f = poly([-1,0,1])
    kappa = 1
    ab = r_jacobi(1+p,0,0)
    xw = gauss(ab)
    X = linspace(a,b,n+1)
    els,G,x,phi = fea_diri2(X, p, kappa, f, xw)

    pp = ppolyfea1sol(els,G,x,phi)
    pf = polytoppoly(f,els,a,b)
    perr = ppolyaxpy(-1,pp,pf)
    ploterr([perr])
    #xx,yy = ppolyval(pp,50)
    #pl.plot(xx,yy)


plotels(X)
#plotbasis(phi)

pl.axis([-1,1,-0.005,0.005])

pl.show()
quit()

ax = pl.subplot(111)
pl.subplots_adjust(bottom=0.25)

pl.axis([-1,1,-2,2])
pl.grid(True)
ps = Slider(pl.axes([0.1, 0.1, 0.7, 0.03]), 'n', 2, 15, valinit=2, valfmt='%1.0f')

def update(val):
    ax.clear()
    p = int(ps.val)
    X = chebyx(p)
    L = lagrange(X)
    xx = pl.linspace(-1, 1, 100)
    ax.plot(X,pl.zeros(p),'o')
    ax.plot(X,pl.ones(p),'o')
    for i in range(p):
        ax.plot(xx,polyvalv(L[i,:],xx))
    pl.draw()
    
ps.on_changed(update)

pl.title('Lobatto basis polynomials')
pl.show()

