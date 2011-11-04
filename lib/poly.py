from matrix_utils import *
from recurrence import r_jacobi
import copy

def poly(L):
    return matrix(L).T

# Convert recurrence coefficients to polynomials
def ab_to_poly(ab):
    n = ab.rows
    P = zeros(n)
    P[0,n-1] = mpf(1)
    P[1,:] = row_join(zeros(1,n-2), matrix([mpf(1), -ab[0,0]]).T)
    for k in range(2,n+1):
        P[k-1,:] = row_join(P[k-2,1:n],matrix([mpf(0)])) - P[k-2,:]*ab[k-2,0] - ab[k-2,1]*P[k-3,:]
    return P

# Convert general recurrence coefficents to polynomials
def r_to_poly(B):
    n = B.rows
    P = zeros(n)
    P[0,n-1] = mpf(1)

    for k in range(n-1):
        S = zeros(1,n)
        for j in range(k+1):
            S = S + B[j,k]*P[k-j,:]
        tpik = row_join(P[k,:], zeros(1))
        P[k+1,:] = tpik[1:n+1] - S
    return P

# Evaluate polynomial on a vector of points (TRACE)
def polyvalv(p,xx):
    yy = zeros(1,len(xx))
    i = 0
    for x in xx:
        yy[i] = polyval(p,x)
        i += 1
    return yy

# Polynomial axpy
def polyaxpy(a,x,y):
    dx,dy = x.cols,y.cols
    if dx < dy:
        z = row_join(zeros(1,dy-dx), a*x)+y
    else:
        z = a*x+row_join(zeros(1,dx-dy), y)
    return z

# Discrete convolution
# From http://www.physics.rutgers.edu/~masud/computing/WPark_recipes_in_python.html
def conv(x, y):
    P, Q, N = len(x), len(y), len(x)+len(y)-1
    z = zeros(1,N)
    for k in range(N):
	t, lower, upper = 0, max(0, k-(Q-1)), min(P-1, k)
	for i in range(lower, upper+1):
	    t = t + x[i] * y[k-i]
	z[k] = t
    return z

# Convolution of multiple polynomials (represented as rows of a matrix P) with a fixed one, p.
def convm(P, p):
    Q = zeros(P.rows, P.cols+len(p)-1)
    for i in range(P.rows):
        Q[i,:] = conv(P[i,:],p)
    return Q

# Convenience routine for monomial x^n
def mon(n):
    p = matrix(1,n+1)
    p[0] = mpf(1)
    return p

# Formal polynomial derivative
def polyder(p):
    N = len(p)
    q = p[0:N-1]
    for i in range(0,N-1): q[i] *= mpf(N-1-i)
    return q

# Formal polynomial primitive function with constant term zero
def polypri(p):
    N = len(p)
    q = row_join(p, zeros(1))
    for i in range(0,N): q[i] /= mpf(N-i)
    return q

# Evaluate the definite integral of the polynomial pq using quadrature rules xw
def quadpq(p, q, xw):
    return fdot(mult(polyvalv(p,xw[:,0]),polyvalv(q,xw[:,0])),xw[:,1])

# Evaluate the definite integral of the polynomial pq on an
# arbitrary interval [a,b] using the quadrature rules xw given on [-1,1]
def quadpqab(p, q, a, b, xw):
    paff,jaff = polyaffj(p,mpf(-1),mpf(1),a,b)
    qaff,jaff = polyaffj(q,mpf(-1),mpf(1),a,b)
    return jaff*fdot(mult(polyvalv(paff,xw[:,0]), polyvalv(qaff,xw[:,0])), xw[:,1])
    
# Compute composition p(s(x))
def polycomp(p, s):
    dp,ds = len(p)-1,len(s)-1 # deg(p),deg(s)
    dps = dp*ds # deg(p(s(x)))=deg(p)deg(s)
    q = zeros(1,dps+1) 
    r = ones(1) 
    for k in range(0,dp+1):
        q = q + row_join(zeros(1,dps-k*ds), p[dp-k]*r) # deg(r)=k*deg(s)
        r = conv(r,s) # r=s^k
    return q

# Compose polynomial p with the affine transform [a,b]->[c,d]
def polyaff(p, a, b, c, d):
    return polycomp(p, matrix([(c-d)/(a-b), (a*d-b*c)/(a-b)]))

## def polyaff0(p, a, b, c, d):
##     aff = matrix([(c-d)/(a-b), (a*d-b*c)/(a-b)])
##     n = len(p)
##     q = row_join(zeros(1,n-1), matrix([p[n-1]]))
##     r = matrix([mpf(1)])
##     for k in range(n-1,0,-1):
##         r = conv(r,aff)
##         q = q + row_join(zeros(1,k-1), p[k-1]*r);
##     return q

# Same as polyaff but also returns the Jacobian of the transform
def polyaffj(p, a, b, c, d):
    return polyaff(p, a, b, c, d), ((c-d)/(a-b))

# affine transform quadrature rule
def xwaff(xw,a,b,c,d):
    xw1 = copy.copy(xw)
    for i in range(xw.rows):
        xw1[i,0] = xw[i,0]*(c-d)/(a-b)+(a*d-b*c)/(a-b)
        xw1[i,1] = xw[i,1]*(c-d)/(a-b)
    return xw1

# Calculate L2 norm squared using quadrature rules xw
def polyl2(p, xw):
    return quadpq(p,p,xw)

# Calculate L2 norm squared using quadrature rules xw on the interval [a,b]
# The quadrature rules are assumed to have support [-1,1]
def polyl2ab(p, a, b, xw):
    paff = polyaff(p,mpf(-1),mpf(1),a,b)
    return quadpq(paff,paff,xw)

# L2 norm squared
def l1norm2(p,xw):
    return quadpq(p,p,xw)

# L2 norm squared on [a,b]
def l1norm2ab(p,a,b,xw):
    return quadapq(p,p,a,b,xw)

# H1 inner product k*\int pq+l*\int p'q'
def h1inner(p,q,xw,k0=mpf(1),k1=mpf(1)):
    i0 = quadpq(p,q,xw)
    i1 = quadpq(polyder(p),polyder(q),xw)
    return k0*i0+k1*i1

# H1 inner product k*\int pq+l*\int p'q' on [a,b]
def h1innerab(p,q,a,b,xw,k0=mpf(1),k1=mpf(1)):
    i0 = quadpqab(p,q,a,b,xw)
    i1 = quadpqab(polyder(p),polyder(q),a,b,xw)
    return k0*i0+k1*i1

# H1 norm squared
def h1norm2(p,xw,k=mpf(1),l=mpf(1)):
    return h1inner(p,p,xw,k,l)

# H1 norm squared on [a,b]
def h1norm2ab(p,a,b,xw,k=mpf(1),l=mpf(1)):
    return h1innerab(p,p,a,b,xw,k,l)

# Chebyshev nodes on [-1,1]
def cheby1(n):
    l = lambda k: cos(pi*mpf(2*k-1)/mpf(2*n))
    return matrix([map(l, range(n,0,-1))])

# Extended Chebyshev nodes of the second kind on [-1,1]
def cheby2(n):
    l = lambda k: cos(pi*mpf(n-k)/mpf(n))
    return matrix([map(l, range(0,n+1))])

# Compute Lagrange basis polynomials on nodes X
def lagrange(X):
    n = len(X)
    if n < 1:
        raise ValueError, "Invalid nodes"

    L = zeros(n,n)
    for j in range(n):
        Lj = [mpf(1)]
        denom = mpf(1)
        for k in range(n):
            if k != j:
                Lj = conv(Lj, [mpf(1), -X[k]])
                denom *= X[j]-X[k]
        L[j,:] = Lj/denom

    return L

def lobatto(n):
    ab = r_jacobi(n,0,0)
    leg = ab_to_poly(ab)
    lob = zeros(n+1)
    lob[0,:] = row_join(zeros(1,n-1), matrix([[mpf(-0.5), mpf(0.5)]]))
    lob[1,:] = row_join(zeros(1,n-1), matrix([[mpf(0.5), mpf(0.5)]]))

    for k in range(2,n):
        leg[k-1,:] /= polyval(leg[k-1,:], 1)
        lob[k,:] = polypri(leg[k-1,:])
        lob[k,n] = -polyval(lob[k,:],-1)
        lob[k,:] *= sqrt(mpf(2*k-1)/mpf(2))
    
    return lob

