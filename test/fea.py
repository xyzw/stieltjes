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

from matplotlib import rc
rc('text', usetex=True)
#rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':'18'})


# Plot element bars
def plotels(X,color="grey",width=1):
    for l in range(len(X)):
        pl.plot([X[l],X[l]],[-1000,1000],color=color,linestyle='--',linewidth=width)
        
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
    dexact = spikeder(s)
    rhs = spikerhs(s,kappa2)
    ad = 0
    bd = 0
    titl = r"$(1-x^2)^{{{0:s}}}$".format(nstr(s))
if args.case == "sinc":
    d = 10
    exact = lambda x : sin(pi*d*x)/(pi*d*x) if x != 0 else 1.0
    dexact = lambda x : (pi*d*x*cos(pi*d*x)-sin(pi*d*x))/(pi*d*x**2) if x != 0 else 0
    rhs = rhsmpmathquad(lambda x : (2*pi*d*x*cos(pi*d*x)-(2-pi**2*d**2*x**2)*sin(pi*d*x))/(pi*d*x**3) + kappa2*sin(pi*d*x)/(pi*d*x))
    ad = exact(a)
    bd = exact(b)
    titl = "sinc"
elif args.case == "arctan":
    r = 100
    exact = arctanjump(r)
    dexact = arctanjumpder(r)
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
    dexact = lambda x : 1.0/(1+x+s)
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
elif args.case == "homo":
    ad = 0
    bd = 1
    kappa = sqrt(kappa2)
    c2 = exp(kappa)/(1.0-exp(4.0*kappa))
    c1 = -c2*exp(2*kappa)
    exact = lambda x : c1*exp(kappa*x)+c2*exp(-kappa*x)
    dexact = lambda x: c1*exp(kappa*x)-c2*exp(-kappa*x)
    rhs = lambda phi,r0,r1,el0,el1 : 0.0

    titl = r'$0$'   
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
        iapnel = 3
        iapp = 3

        #xx,yy = ppolyvalres(ppolyder(uh),1000)
        #pl.plot(xx,yy,label=r"$u_h'$")

        duh = ppolyder(uh)
        dduh = ppolyder(duh)

        #xx,yy = ppolyerrvalres(duh,dexact,1000)
        #pl.plot(xx,yy,label=r"$e_h'$")

        #xx,yy = ppolyvalres(dduh,1000)
        #pl.plot(xx,yy,label=r"$u_h''$")

        #res = ppolyaxpy(-kappa2,uh,dduh) # u_h''-kappa2 u_h
        #xx,yy = ppolyaxpbyvalres(1.0,exact,-1.0,res,100)
        #pl.plot(xx,yy,label=r"$\widetilde{{f}}$")

        duX = map(dexact, X)       
        #errh = fea1diap(els,G,x,iapnel,p,iapp,phi,kappa2,0,[ad,bd],duX,rhs)
        #errhhneu,errhhneunorm = iapneu(els,G,x,phi,kappa2,rhs,duX,iapnel,iapp)
        #errhhbub,errhhbubnorm,dt = iapbub(els,G,x,phi,kappa2,rhs,iapnel,iapp)
        #errhhbubxn,errhhbubxnnorm,dt = iapbub(els,G,x,phi,kappa2,rhs,iapnel,iapp,False)
        errhhbubxnspread,errhhbubxnspreadnorm,tbubxnspread = iapbubspread(els,G,x,phi,kappa2,rhs,iapnel,iapp,False)

        #xx,yy = ppolyerrvalres(uh,exact,1000)
        #pl.plot(xx,yy,label=r"$e_h$",linewidth=1.5)

        #xx,yy = ppolyvalres(errhhneu,1000)
        #pl.plot(xx,yy,label=r"$e_{{hh}}^N$",linewidth=2)

        #xx,yy = ppolyvalres(errhhbub,1000)
        #pl.plot(xx,yy,label=r"$e_{{hh}}^B$",linewidth=2)

        #xx,yy = ppolyvalres(errhhbubxn,1000)
        #pl.plot(xx,yy,label=r"$e_{{hh}}^{BM}$",linewidth=2)

        #xx,yy = ppolyvalres(errhhbubxnspread,1000)
        #pl.plot(xx,yy,label=r"$e_{{hh}}^{BMS}$",linewidth=2)

        # Plot norm ratio bars
        #exactnorm = ppolyerrh1norm2intv(uh,exact,dexact)
        #errhhneunorms = ppolyh1norm2intv(errhhneu, gauss(r_jacobi(iapp+1)))
        #for e in range(len(els)):
        #    w = 0.05*(els[e][1]-els[e][0])
        #    pl.bar(0.5*(els[e][0]+els[e][1]-2*w), errhhneunorm[e]/exactnorm[e], w)
        #    pl.bar(0.5*(els[e][0]+els[e][1]+2*w), errhhbubnorm[e]/exactnorm[e], w)
        
        #uh = ppolyrefuni(uh,iapnel)
        #errerrh1 = ppolyerrh1norm2(ppolyaxpy(1,uh,errh),exact,dexact)
        
        #xx,yy = ppolyerrvalres(ppolyaxpy(1,uh,errh),exact,1000)
        #pl.plot(xx,yy,label=r"$e_{{hh}}-e_h$, $\|e_{{hh}}-e_h\|_{{H^1}}^2={0:s}$".format(nstr(errerrh1)))

        for el in els:
            plotels(linspace(el[0],el[1],iapnel+1),"green",1.0)


pl.legend()
pl.title(titl)

plotels(X,"black",2)
#print phi
#plotbasis(phi)

pl.axis([-1,1,ymin,ymax])

pl.show()

