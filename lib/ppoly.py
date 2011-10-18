
from matrix_utils import *
from poly import polyaxpy,polyvalv,polyaff,quadpq,polyl2

#TODO: Improve this

# Piecewise polynomial class
class ppoly(object):
    def __init__(self, I, P=[]):
        self._poly = P if len(P)>0 else [[0] for i in range(len(I))]
        self._intv = I

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
            s += str(self._intv[i]) + " , " + str(self._poly[i].T) + "\n"
        return s

def ppolysamesupp(pp,pq):
     if len(pp.intv) != len(pq.intv): return False
     return True

# Evaluate on the supporting intervals with ires resoltion (i.e. for plotting)
def ppolyval(pp,ires):
    xx = []
    yy = []
    for i in range(len(pp.intv)):
        xxx = linspace(pp.intv[i][0], pp.intv[i][1], ires)
        yyy = polyvalv(pp.poly[i], xxx)
        xx.extend(xxx)
        yy.extend(yyy)
    return xx,yy

## ppoly axpy operation
def ppolyaxpy(a,x,y):
    if not ppolysamesupp(x,y):
        raise ValueError, "ppoly incompatible for axpy"
    z = ppoly(x.intv)
    for i in range(len(z.intv)):
        z.poly[i] = polyaxpy(a, x.poly[i], y.poly[i])
    return z

# Given a poly on intervals I, hull(I)=[a,b] convert it to a ppoly
def polytoppoly(p,I,a,b):
    pp = ppoly(I)
    for i in range(len(I)):
        pp.poly[i] = p
    return pp

# Calculate L2 norm squared of the ppoly, the xw quadrature rules are supported on [-1,1]
def ppolyl2norm(pp, xw):
    norm2 = 0
    for i in range(len(pp.intv)):
        norm2 += polyl2(pp.poly[i], pp.intv[i][0], pp.intv[i][1], xw)
    return norm2

    
