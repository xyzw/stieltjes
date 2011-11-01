from decomp import *

#
# GAUSS
#
# SYNOPSIS
# Generate Gauss quadrature rules from the recursion coefficents of an orthogonal
# polynomial sequence.
#
# INPUT
# ab - Recursion coefficients in the Stieltjes relations as a matrix of size nx2
#
# RETURN VALUE
# Quadrature nodes and weights as an nx2 matrix (usually denoted as xw as per W. Gautschi)
#
# COMMENTS
# This is a Python port of the MATLAB OPQ routine GAUSS. The only difference is
# the use of the more efficient tridiagonal QL algorithm.
#
def gauss(ab):
    if ab.cols != 2:
        raise ValueError, "Not a recurrence coefficient matrix specified"

    n = ab.rows
    d,e = jacobi_matrix(ab)
    V = tql2(d,e,ev=True)
    return row_join(d,ab[0,1]*pwr(V[0,:],2).T)
    
#
# FEJER
#
# SYNOPSIS
# Fejer quadrature rule
#
# INPUT
# N - Point count
#
# RETURN VALUE
# Quadrature nodes and weights as an nx2 matrix (usually denoted as xw as per W. Gautschi)
#
# COMMENTS
# This is a Python port of the MATLAB OPQ routine FEJER. 
#
def fejer(N):
    n = matrix(range(N,0,-1))
    th = matrix(map(lambda k: mpf(2*k-1)*pi/mpf(2*N), n))
    uv = zeros(N,2)
    uv[:,0] = matrix(map(lambda x: cos(x), th))
    for k in range(N,0,-1):
        s = sum([cos(2*m*th[k-1])/mpf(4*(m**2)-1) for m in range(1,int(floor(N/2))+1)])
        uv[k-1,1]=mpf(2*(1-2*s))/mpf(N)
    return uv

#
# MCDIS
#
# SYNOPSIS
# Multiple component measure discretization
#
# INPUT
# N - Point count
#
# RETURN VALUE
# Quadrature nodes and weights as an nx2 matrix (usually denoted as xw as per W. Gautschi)
#
# COMMENTS
# This is a Python port of the MATLAB OPQ routine MCDIS. 
#
def mcdis(N,eps0,quad,Nmax):
   # global mc mp iq idelta irout DM uv AB
    if N<1:
        print "Input variable N out of range"
        return
    mc = 5
    idelta=1
    iNcap=1
    kount=-1
    ab = zeros(N,2)
    b = ones(N,1)
    Ncap = floor((2*N-1)/idelta)
    while any(fabs(ab[k,1]-b[k]) > eps0*fabs(ab[k,1]) for k in range(N)):
          b = ab[:,2]
          kount = kount+1
          if kount > 1: iNcap=2**(int(floor(kount/5)))*N
          Ncap = Ncap + iNcap
          if Ncap > Nmax:
              print "Ncap exceeds Nmax in mcdis with irout=", irout
              return
          
          mtNcap = mc*Ncap
          if quad is None: uv = fejer(Ncap)
          for i in range(1,mc+1):
              im1tn = (i-1)*Ncap
              xw = quad(Ncap,i)
              xwm[im1tn+1:im1tn+Ncap,:] = xw

          ab = stieltjes(N,xwm)
          
    return ab
