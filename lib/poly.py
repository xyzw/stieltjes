from matrix_utils import *

# Convert recurrence coefficients to polynomials
def ab_to_poly(ab):
    n = ab.rows
    P = zeros(n)
    P[0,n-1] = mpf(1)
    P[1,:] = row_join(zeros(1,n-2), matrix([mpf(1), -ab[0,0]]).T)
    for k in range(3,n+1):
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
