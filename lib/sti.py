from matrix_utils import *

#
# DSTI
#
# SYNOPSIS
# Discrete Stieltjes algorithm for the generation of recursion coeffcients
# from quadrature rules.
#
# INPUT
# xw - An nx2 matrix containing nodes and weights.
#
# RETURN VALUE
# Recursion coefficents in the form of an nx2 matrix.
#
# COMMENTS
# This is a Python port of the OPQ MATLAB routine STIELTJES.
#
def dsti(xw):
    if xw.cols != 2:
        raise ValueError, "Not a quadrature rule specified"
    
    n = xw.rows
    x,w = xw[:,0], xw[:,1]
    ab = zeros(n,2)
    s0 = fsum(w)
    ab[0,0] = fdot(x,w)/s0
    ab[0,1] = s0

    if n == 1: return ab

    p1 = zeros(n,1)
    p2 = ones(n,1)

    for k in range(n-1):
        p0 = p1
        p1 = p2
        p2 = mult(x-ab[k,0],p1)-ab[k,1]*p0

        s1 = fdot(w,pwr(p2,2))
        s2 = fdot(x,mult(w,pwr(p2,2)))
        
        ab[k+1,0],ab[k+1,1] = s2/s1, s1/s0
        s0 = s1

    return ab
    
#
# DSTIS1
#
# SYNOPSIS
# Discrete Stieltjes algorithm for the generation of recursion coeffcients
# of a first-order Sobolev orthogonal sequence from the quadrature rules
# of the Sobolev-inner product.
#
# INPUT
# xw - An nx2 matrix containing nodes and weights.
#
# RETURN VALUE
# Recursion coefficents in the form of an nx2 matrix.
#
# COMMENTS
# This is a Python port of the OPQ MATLAB routine STIELTJES_SOB,
# with the simplification that s=1. 
#
def dstis1(xw,nd,a0):
    if xw.cols != 4:
        raise ValueError, "Must specify two quadrature rules in parameter xw"

    n = xw.rows
    md = max(nd)
    B = zeros(n)

    # Initialization
    p = [zeros(n,md), zeros(n,md)]
    p[0][0,:] = mpf(1)
    xp = [[zeros(n,md), zeros(n,md)], [zeros(n,md), zeros(n,md)]]

    for s in [0,1]:
        for ss in [0,1]:
            for mu in range(nd[ss]):
                if s == 0:
                    xp[s][ss][0,mu] = xw[mu,ss]
                elif s == 1:
                    xp[s][ss][0,mu] = mpf(1)
                else:
                    raise AssertError, "WHF?"

    B[0,0] = a0
    if n == 1: return

    # Continuation
    for k in range(2,n+1):
        for j in range(1,k+1):
            for s in [0,1]:
                for nu in range(nd[s]):
                    sm = xp[s][s][k-2,nu]
                    for l in range(1,k):
                        sm = sm - B[l-1,k-2] * p[s][k-l-1,nu]
                    p[s][k-1,nu] = sm


            Bnum = mpf(0)
            Bden = mpf(0)
            
            for s in [0,1]:
                for ss in [0,1]:
                    for mu in range(nd[ss]):
                        xsm = xw[mu,ss]*xp[s][ss][k-2,mu]

                        if s == 1: xsm = xsm + xp[0][ss][k-2,mu]

                        for l in range(1,k):
                            xsm = xsm - B[l-1,k-2]*xp[s][ss][k-l-1,mu]

                        xp[s][ss][k-1,mu] = xsm

                        if s == ss:
                            Bnum = Bnum + xw[mu,s+2]*xp[s][ss][k-1,mu]*p[s][k-j,mu]
                            Bden = Bden + xw[mu,s+2]*(p[s][k-j,mu]*p[s][k-j,mu])

            B[j-1,k-1] = Bnum/Bden
            #return B
                        
    return B
    
#
# PDSTIS1
#
# SYNOPSIS
# Discrete Stieltjes algorithm for the generation of recursion coeffcients
# of a first-order Sobolev orthogonal sequence from the quadrature rules
# of the Sobolev-inner product.
#
# INPUT
# xw - An nx2 matrix containing nodes and weights.
#
# RETURN VALUE
# Recursion coefficents in the form of an nx2 matrix.
#
# COMMENTS
# This is a Python port of the OPQ MATLAB routine STIELTJES_SOB,
# with the simplification that s=1. 
#
def pdstis1(xw,nd,a0):
    if xw.cols != 4:
        raise ValueError, "Must specify two quadrature rules in parameter xw"

    print xw

    n = xw.rows
    md = max(nd)
    B = zeros(n)

    # Initialization
    p = [zeros(n,md), zeros(n,md)]
    p[0][0,:] = mpf(1)
    xp = [[zeros(n,md), zeros(n,md)], [zeros(n,md), zeros(n,md)]]

    for s in [0,1]:
        for ss in [0,1]:
            for mu in range(nd[ss]):
                if s == 0:
                    xp[s][ss][0,mu] = xw[mu,ss]
                elif s == 1:
                    xp[s][ss][0,mu] = mpf(1)
                else:
                    raise AssertError, "WHF?"

    B[0,0] = a0
    if n == 1: return

    # Continuation
    for k in range(2,n+1):
        for j in range(1,k+1):
            for s in [0,1]:
                for nu in range(nd[s]):
                    sm = xp[s][s][k-2,nu]
                    for l in range(1,k):
                        sm = sm - B[l-1,k-2] * p[s][k-l-1,nu]
                    p[s][k-1,nu] = sm


            Bnum = mpf(0)
            Bden = mpf(0)
            
            for s in [0,1]:
                for ss in [0,1]:
                    for mu in range(nd[ss]):
                        xsm = xw[mu,ss]*xp[s][ss][k-2,mu]

                        if s == 1: xsm = xsm + xp[0][ss][k-2,mu]

                        for l in range(1,k):
                            xsm = xsm - B[l-1,k-2]*xp[s][ss][k-l-1,mu]

                        xp[s][ss][k-1,mu] = xsm

                        if s == ss:
                            Bnum = Bnum + xw[mu,s+2]*xp[s][ss][k-1,mu]*p[s][k-j,mu]
                            Bden = Bden + xw[mu,s+2]*(p[s][k-j,mu]*p[s][k-j,mu])

            B[j-1,k-1] = Bnum/Bden
            #return B
                        
    return B
    
