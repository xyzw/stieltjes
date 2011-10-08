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
# Quadrature nodes and weights as an nx3 matrix (usually denoted as xw as per W. Gautschi)
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
    
