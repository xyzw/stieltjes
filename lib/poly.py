from matrix_utils import *
from recurrence import r_jacobi

def poly(L):
    return matrix(L)

# Convert recurrence coefficients to polynomials
def ab_to_poly(ab):
    n = ab.rows
    P = zeros(n)
    P[0,n-1] = mpf(1)
    P[1,:] = row_join(zeros(1,n-2), matrix([mpf(1), -ab[0,0]]).T)
    for k in range(2,n+1):
        P[k-1,:] = row_join(P[k-2,1:n],matrix([mpf(0)])) - P[k-2,:]*ab[k-2,0] - ab[k-2,1]*P[k-3,:]
        
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
    dx,dy = len(x),len(y)
    if dx < dy:
        z = col_join(zeros(1,dy-dx), a*x)+y
    else:
        z = a*x+col_join(zeros(1,dx-dy), y)
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

# Formal polynomial derivative
def polyder(p):
    N = len(p)
    q = p[0:N-1]
    for i in range(0,N-1): q[i]*= N-1-i
    return q

# Formal polynomial primitive function with constant term zero
def polypri(p):
    N = len(p)
    q = row_join(p, zeros(1))
    for i in range(0,N): q[i] /= N-i
    return q

# Evaluate integral of the polynomial pq using quadrature rules xw
def quadpq(p, q, xw):
    return fdot(mult(polyvalv(p,xw[:,0]),polyvalv(q,xw[:,0])),xw[:,1])

# Compose polynomial p with the affine transform [a,b]->[c,d]
def polyaff(p, a, b, c, d):
    aff = matrix([(c-d)/(a-b), (a*d-b*c)/(a-b)])
    n = len(p);
    q = row_join(zeros(1,n-1), matrix([p[n-1]]))
    r = matrix([mpf(1)])
    for k in range(n-1,0,-1):
        r = conv(r,aff)
        q = q + row_join(zeros(1,k-1), p[k-1]*r);
    return q

# Chebyshev nodes on [-1,1]
def cheby(n):
    l = lambda k: cos(pi*mpf(2*k-1)/mpf(2*n))
    return matrix([map(l, range(n,0,-1))])

# Extended Chebyshev nodes on [-1,1]
def chebyx(n):
    l = lambda k: cos(pi*mpf(2*k-1)/mpf(2*n)) / cos(pi/mpf(2*n))
    return matrix([map(l, range(n,0,-1))]) # reversed

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
                Lj = conv(Lj, [mpf(1), X[k]])
                denom *= X[k]-X[j]
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
