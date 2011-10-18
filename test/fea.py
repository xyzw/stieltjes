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
from matplotlib.pyplot import legend

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


pl.grid(True)

mp.dps = 40

kappa = 1
a=-1
b=1
n=5

s = 25
spike = lambda x : (1-x**2)**s
spikerhs = lambda x : 2*s*(1-x**2)**(s-1)-s*(s-1)*4*x**2*(1-x**2)**(s-2) + (1-x**2)**s

for p in range(1,3):


    ab = r_jacobi(1+p,0,0)
    xw = gauss(ab)
    X = linspace(a,b,n+1)

    rhs = rhsmpmathquad(spikerhs)
    els,G,x,phi = fea_diri2(X, p, kappa, [0,0], rhs, xw)

    pp = ppolyfea1sol(els,G,x,phi)
    pf = polytoppoly(poly([-1,0,1]),els,a,b)
    perr = ppolyaxpy(-1,pp,pf)
    #ploterr(p,[perr])
    #print "l2err =", ppolyl2norm(perr,xw)

    xx,yy = ppolyval(pp,50)
    pl.plot(xx,yy,label="p=" + str(p))


xx = linspace(a,b,1000)
yy = map(spike, xx)
pl.plot(xx,yy,linewidth=2,label="exact")

pl.title("An equation of which the exact solution is (1-x**2)**s (n=" + str(n) + ", s=" + str(s) + ")")
pl.legend()

plotels(X)
#print phi
#plotbasis(phi)

pl.axis([-1,1,-2,2])

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

