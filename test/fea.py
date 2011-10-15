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

x = matrix([1,1,1])
y = matrix([1,1,1,1])
print x.T, y.T
print polyaxpy(-10,x,y).T

pp = ppoly([[0,1],[1,2],[2,3]], [poly([1]), poly([2]), poly([3])])
pq = ppoly([[0,1],[1,2],[2,3]], [poly([1,0]), poly([2,0]), poly([3,0])])

print pp
print pq

print ppolyaxpy(1,pp,pq)

quit()

mp.dps = 10
n=4
p=6
kappa = 1
ab = r_jacobi(n,0,0)
xw = gauss(ab)

X = linspace(-1,1,n+1)
els,G,x,phi = fea_diri2(X, p, kappa, matrix([-1,0,1]), xw)

#print ab_to_poly(ab)

f = matrix([-1,0,1])
#ftilde = polyaff(f, -0.5, 0.5, -1, 1)
#xx = pl.arange(-1.0, 1.0, 0.01)

#pl.plot(xx,polyvalv(ftilde,xx))
#for i in range(p+1):
 #   pl.plot(xx,polyvalv(phi[i,:],xx))
xx,yy = evalfea1sol(els,G,x,phi,100)
pl.plot(xx,yy)
pl.plot(xx,polyvalv(f,xx))
for l in range(len(X)):
    pl.plot([X[l],X[l]],[-1000,1000],color='grey',linestyle='--')
pl.axis([-1,1,-2,2])
pl.grid(True)
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

