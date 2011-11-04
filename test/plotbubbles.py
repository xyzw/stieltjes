import sys
sys.path.append('../lib')
import time
import argparse

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
#rc('text', usetex=True)
#rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':'18'})


# Plot element bars
def plotels(X):
    for l in range(len(X)):
        pl.plot([X[l],X[l]],[-1000,1000],color='grey',linestyle='--')
        

mp.dps = 20

parser = argparse.ArgumentParser(description='Plot bubble functions.',add_help=False)
parser.add_argument('-n', metavar='NUMBUB', type=int, required=True, help='Number of bubbles.')
parser.add_argument('-h', metavar='INTV', type=float, required=True, help='Size of support interval.')
parser.add_argument('-kappa2', metavar='KAPPA2', type=float, required=True, help='Constant in DE.')
parser.add_argument('-alpha', metavar='ALPHA', type=int, required=True, help='Weight function degree.')
parser.add_argument('-p', metavar='FEDEG', type=int, required=True, help='Finite element degree.')
parser.add_argument('-nel', metavar='NEL', type=int, required=True, help='Finite elements.')
args = parser.parse_args()

n = args.n
h = mpf(args.h)
kappa2 = mpf(args.kappa2)
alpha = args.alpha
p = args.p
nel = args.nel


Chat = bubble(n,h,kappa2,alpha,p,nel)

for k in range(n):
    xx,yy = ppolyval(Chat[k],50)
    pl.plot(xx,yy,label=r"$C_{0:d}$".format(k),linewidth=0.8)

pl.title(r"$\{{\widehat{{C}}_k\}}$ on $[{4:s},{5:s}]$ ($\kappa^2={0:s}$, $n={1:d}$, $\alpha={2:d}$, $p={3:d}$)".format(nstr(kappa2),n,alpha,p,nstr(-mpf(0.5)*h),nstr(mpf(0.5)*h)))
pl.legend()

plotels(linspace(-mpf(0.5)*h,mpf(0.5)*h,nel+1))

pl.axis([float(-0.5*h),float(0.5*h),-0.2,0.2])
pl.grid(True)

pl.show()
