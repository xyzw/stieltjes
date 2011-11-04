import sys
sys.path.append('../lib')
import argparse
from recurrence import *
from decomp import *
from quad import *
from sti import *
from itertools import product
from poly import *
from ppoly import *
from fea1d import *
import pylab as pl
from matplotlib.widgets import Slider
from matplotlib.pyplot import legend
from egrhs import *
from aposteriori import *


# Plot element bars
def plotels(X,color="grey"):
    for l in range(len(X)):
        pl.plot([X[l],X[l]],[-1000,1000],color=color,linestyle='--')
        
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

parser = argparse.ArgumentParser(description='Plot finite element problems.',add_help=False)
parser.add_argument('-kappa2', metavar='KAPPA2', type=float, required=True, help='Constant in DE.')
parser.add_argument('-nel', metavar='NEL', type=int, required=True, help='Finite elements.')
parser.add_argument('-pmin', metavar='PMIN', type=int, required=True, help='Finite element minimum degree')
parser.add_argument('-pmax', metavar='PMAX', type=int, required=True, help='Finite element maximum degree')
parser.add_argument('-case', metavar='PARAM', required=True, help='Test case name.')
parser.add_argument('-mode', metavar='MODE', required=True, help='solve, iap')
parser.add_argument('-ymin', metavar='YMIN', type=float, required=True, help='Plot ordinate interval left endpoint.')
parser.add_argument('-ymax', metavar='YMAX', type=float, required=True, help='Plot ordinate interval right endpoint.')
args = parser.parse_args()

kappa2 = args.kappa2
a=-1
b=1
nel=args.nel
pmin=args.pmin
pmax=args.pmax
ymin=args.ymin
ymax=args.ymax

titl = None
bt = 0

if args.case == "spike":
    s = 100
    exact = spike(s)
    rhs = spikerhs(s,kappa2)
    ad = 0
    bd = 0
    titl = r"$(1-x^2)^{{{0:s}}}$".format(nstr(s))
elif args.case == "arctan":
    r = 500
    exact = arctanjump(r)
    rhs = arctanjumprhs(r,kappa2)
    ad = exact(a)
    bd = exact(b)
    titl = r'$\frac{{2}}{{\pi}}\arctan({0:s}x)$'.format(nstr(r))
elif args.case == "arctanbub":
    r = 100
    exact = arctanbubjump(r)
    rhs = arctanbubjumprhs(r,kappa2)
    ad = exact(a)
    bd = exact(b)
    titl = r'$\frac{{2}}{{\pi}}\arctan({0:s}x)(1-x^2)$'.format(nstr(r))
elif args.case == "expspike":
    s = 0.01
    exact = expspike(s)
    rhs = expspikerhs(s,kappa2)
    ad = exact(a)
    bd = exact(b)
    titl = r"$e^{{-x^2/{0:s}^2}}$".format(nstr(s))
elif args.case == "loggap":
    s = 0.001
    exact = loggap(s)
    rhs = loggaprhs(s,kappa2)
    ad = exact(a)
    bd = exact(b)
    titl = r"$\log(1+x+{0:s})$".format(nstr(s))
elif args.case == "neuspike":
    s = 100
    bt = 1
    exact = neuspike(s)
    rhs = neuspikerhs(s,kappa2)
    ad = neuspikeder(s)(a)
    bd = neuspikeder(s)(b)
    titl = r"$(1-x^2)^{{{0:s}}}$".format(nstr(s))
elif args.case == "neuarctan":
    r = 500
    bt = 1
    exact = arctanjump(r)
    rhs = arctanjumprhs(r,kappa2)
    ad = arctanjumpder(r)(a)
    bd = arctanjumpder(r)(b)
    titl = r'$\frac{{2}}{{\pi}}\arctan({0:s}x)$'.format(nstr(r))
elif args.case == "neuarctanbub":
    r = 100
    bt = 1
    exact = arctanbubjump(r)
    rhs = arctanbubjumprhs(r,kappa2)
    ad = arctanbubjumpder(r)(a)
    bd = arctanbubjumpder(r)(b)
    titl = r'$\frac{{2}}{{\pi}}\arctan({0:s}x)(1-x^2)$'.format(nstr(r))
else:
    print "error: invalid test case specified. specify one of the following."
    print "spike arctan arctanbub expspike loggap neuspike"
    quit()

X = linspace(a,b,nel+1)

for p in range(pmin,pmax+1):
    ab = r_jacobi(p+1)
    xw = gauss(ab)

    print ">>>> Finite element analysis for p=", p
    phi = lagrangecheby(p)
    els,G,x,phi = fea1dh(X, lagrangecheby(p), kappa2, bt, [ad,bd], rhs)
    uh = ppolyfea1sol(els,G,x,phi)

    if args.mode == "solve":
        xx,yy = ppolyvalres(uh,50)
        pl.plot(xx,yy,label=r"$u_h$")

        duh = ppolyder(uh)
        xx,yy = ppolyvalres(duh,1000)
        pl.plot(xx,yy,label=r"$u_h'$")
        
        xx = linspace(a,b,1000)
        yy = map(exact, xx)
        pl.plot(xx,yy,linewidth=2,label=r"$u$")
    elif args.mode == "iap":
        print ">>>> Implicit a posteriori error estimation"
        iapnel = 2
        iapp = 5
        
        errh = fea1diap(els,G,x,iapnel,iapp,phi,kappa2,0,[ad,bd],rhs)
        
        xx,yy = ppolyerrvalres(uh,exact,1000)
        pl.plot(xx,yy,label=r"$e$")

        errhh1 = ppolyh1norm2(errh, gauss(r_jacobi(iapp+1)))
        xx,yy = ppolyvalres(errh,1000)
        pl.plot(xx,yy,label=r"$e_h$, $\|e_h\|_{{H^1}}={0:s}$".format(nstr(errhh1)))
        
        #xx,yy = ppolyvalres(ppolyder(uh),1000)
        #pl.plot(xx,yy,label=r"$u_h'$")
        
        for el in els:
            plotels(linspace(el[0],el[1],iapnel+1),"green")


pl.legend()
pl.title(titl)

plotels(X,"black")
#print phi
#plotbasis(phi)

pl.axis([-1,1,ymin,ymax])

pl.show()

