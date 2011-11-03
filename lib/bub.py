from poly import *
from chebyshev import *
from sti import *
from quad import *

# Bubble weight (1-x^(2*alpha))
def bw(alpha=1):
    w = -mon(2*alpha)
    w[2*alpha] = mpf(1)
    return w

# Calculate modified moments on [a,b]
def bwmom(n,a,b,kappa2,alpha,xwl):
    w = polyaff(bw(alpha),a,b,mpf(-1),mpf(1))
    mu0 = polyaxpy(kappa2,conv(w,w),conv(polyder(polyder(w)),-w)) # k^2w^2-w''w
    mu1 = conv(w,w)
    mom0 = zeros(1,2*n);
    mom1 = zeros(1,2*n);
    for k in range(0,2*n):
        mom0[k] = quadpqab(mon(k),mu0,a,b,xwl) # int x^k dmu0
        mom1[k] = quadpqab(mon(k),mu1,a,b,xwl) # int x^k dmu1
    return mom0,mom1

def prebub(n,h,kappa2,alpha,xwl):
    a = -mpf(0.5)*h
    b = mpf(0.5)*h

    mom0,mom1 = bwmom(n,a,b,kappa2,alpha,xwl)
    ab0,n2 = chebyshev(n,mom0)
    ab1,n2 = chebyshev(n,mom1)
    xw0 = gauss(ab0)
    xw1 = gauss(ab1)
    XW = row_join( row_join(xw0[:,0],xw1[:,0]), row_join(xw0[:,1],xw1[:,1]) )
    B,normsq = dstis1(n,XW,[n, n],0)

    P = r_to_poly(B)
    bwaff = polyaff(bw(alpha),a,b,mpf(-1),mpf(1))
    P = convm(P, bwaff)

    # H1-normalize
    for k in range(n):
        h1ni = mpf(1)/sqrt(h1norm2ab(P[k,:],a,b,xwl))
        P[k,:] *= h1ni

    return P

def gramh1(P,a,b,kappa2,xwl):
    n = P.rows
    G = zeros(n)
    for k in range(n):
        for l in range(n):
            G[k,l] = h1innerab(P[k,:],P[l,:],a,b,xwl,kappa2)
    return G

# Compute the absolute maximum of off-diagonal elements
def maxod(G):
    n = G.rows
    maxi = fabs(G[0,1])
    for k in range(n):
        for l in range(n):
            if k != l:
                maxi = max(maxi, fabs(G[k,l]))
    return maxi
