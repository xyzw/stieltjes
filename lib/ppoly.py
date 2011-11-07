
from matrix_utils import *
from poly import *
import copy

#TODO: Improve this

#
# ppoly - piecewise polynomial class
#
class ppoly(object):
    def __init__(self, I=[], P=[]):
        self._poly = P if len(P)>0 else [zeros(1) for i in range(len(I))]
        self._intv = copy.copy(I)

    @property
    def poly(self):
        return self._poly

    @poly.setter
    def poly(self, p):
        self._poly = p

    @property
    def intv(self):
        return self._intv

    @intv.setter
    def intv(self,i):
        self._intv = i

    def __str__(self):
        s = ""
        for i in range(len(self._intv)):
            s += str(self._intv[i]) + " , " + str(self._poly[i]) + "\n"
        return s

def ppolysamesupp(pp,pq):
     if len(pp.intv) != len(pq.intv): return False
     return True

# Evaluate on the supporting intervals with ires resoltion (i.e. for plotting)
def ppolyvalres(pp,ires):
    xx = []
    yy = []
    for i in range(len(pp.intv)):
        xxx = linspace(pp.intv[i][0], pp.intv[i][1], ires)
        yyy = polyvalv(pp.poly[i], xxx)
        xx.extend(xxx)
        yy.extend(yyy)
    return xx,yy

def ppolyval(pp,x):
    for i in range(len(pp.intv)):
        if pp.intv[i][0] <= x and x <= pp.intv[i][1]:
            return polyval(pp.poly[i],x)

    raise ValueError, "Specified point does not lie in the support"

def ppolyzero(intv):
    pp = ppoly(intv)
    for i in range(len(intv)):
        pp.poly[i] = zeros(1)
    return pp

# Extend pp with qq
def ppolyext(pp, qq):
    pp.intv.extend(qq.intv)
    pp.poly.extend(qq.poly)

## ppoly axpy operation
def ppolyaxpy(a,x,y):
    if not ppolysamesupp(x,y):
        raise ValueError, "ppoly incompatible for axpy"
    z = ppoly(y.intv)
    for i in range(len(z.intv)):
        z.poly[i] = polyaxpy(a, x.poly[i], y.poly[i])
    return z

# Refines pp's intvervals into n subintervals uniformly
def ppolyrefuni(pp,n):
    rr = ppoly()
    for i in range(len(pp.intv)):
        X = linspace(pp.intv[i][0],pp.intv[i][1],n+1)
        refintv = zip(X[0:n], X[1:n+1])
        rr.intv.extend(refintv)
        for intv in refintv:
            rr.poly.append(pp.poly[i])
    return rr

# Given a poly on intervals I convert it to a ppoly
def polytoppoly(p,I):
    pp = ppoly(I)
    for i in range(len(I)):
        pp.poly[i] = p
    return pp

def ppolyder(pp):
    dpp = ppoly(pp.intv)
    for i in range(len(pp.poly)):
        dpp.poly[i] = polyder(pp.poly[i])
    return dpp

def ppolyaff(pp,a,b,c,d):
    alpha = (c-d)/(a-b)
    beta = (a*d-b*c)/(a-b)
    qq = ppoly(pp.intv)
    for i in range(len(pp.poly)):
        qq.poly[i] = polyaff(pp.poly[i],a,b,c,d)
        qq.intv[i] = (alpha*pp.intv[i][0] + beta, alpha*pp.intv[i][1] + beta)
    return qq

def ppolyscale(pp,a):
    qq = ppoly(pp.intv)
    for i in range(len(pp.intv)):
        qq.poly[i] = a*pp.poly[i]
    return qq


def quadppqq(pp,qq,xw):
    s = 0
    for i in range(len(pp.intv)):
        if pp.intv[i][0] != pp.intv[i][0] or pp.intv[i][1] != pp.intv[i][1]:
            raise ValueError, "Interval mismatch"
        
        ppaff = polyaff(pp.poly[i], mpf(-1), mpf(1), pp.intv[i][0], pp.intv[i][1])
        qqaff = polyaff(qq.poly[i], mpf(-1), mpf(1), qq.intv[i][0], pp.intv[i][1])
        s += quadpq(ppaff, qqaff, xw)

    return s

# Compute the difference of a ppoly and a general function f and evaluate it just like ppolyvalres
def ppolyerrvalres(pp, f, ires):
    xx = []
    yy = []
    for i in range(len(pp.intv)):
        xxx = linspace(pp.intv[i][0], pp.intv[i][1], ires)
        yyy = matrix([[f(x) for x in xxx]])-polyvalv(pp.poly[i], xxx)
        xx.extend(xxx)
        yy.extend(yyy)
    return xx,yy

def ppolyaxpbyvalres(a, f, b, pp, ires):
    xx = []
    yy = []
    for i in range(len(pp.intv)):
        xxx = linspace(pp.intv[i][0], pp.intv[i][1], ires)
        yyy = a*matrix([[f(x) for x in xxx]]) + b*polyvalv(pp.poly[i], xxx)
        xx.extend(xxx)
        yy.extend(yyy)
    return xx,yy


# Compute H1 norm squared of an error
def ppolyerrh1norm2intv(pp, f, df, k=mpf(1), l=mpf(1)):
    s = []
    for i in range(len(pp.intv)):
        integ = lambda x : k*(polyval(pp.poly[i],x) - f(x))**2 + l*(polyval(polyder(pp.poly[i]),x) - df(x) )**2
        s.append(quad(integ, pp.intv[i]))
    return s

# Compute H1 inner product
def ppolyh1inner(pp,qq,xw,k=mpf(1),l=mpf(1)):
    s = 0
    for i in range(len(pp.intv)):
        s += h1innerab(pp.poly[i], qq.poly[i], pp.intv[i][0], pp.intv[i][1], xw, k, l)
    return s

# Compute H1 norm squared
def ppolyh1norm2(pp,xw,k=mpf(1),l=mpf(1)):
    s = 0
    for i in range(len(pp.intv)):
        s += h1norm2ab(pp.poly[i], pp.intv[i][0], pp.intv[i][1], xw, k, l)
    return s

def ppolyh1norm2intv(pp,xw,k=mpf(1),l=mpf(1)):
    norms = []
    for i in range(len(pp.intv)):
        norms.append(h1norm2ab(pp.poly[i], pp.intv[i][0], pp.intv[i][1], xw, k, l))
    return norms
    
