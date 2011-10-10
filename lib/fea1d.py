from matrix_utils import *
from poly import *
from itertools import product


# Subdivide [a,b] to n uniform elements
def equidistant(a,b,n):
    if b<=a:
        raise ValueError, "Invalid interval specfied"

    if n<=0:
        raise ValueError, "Nonpositive partition size requested"
    
    h = mpf(b-a)/mpf(n)
    x0 = arange(a,b-h,h)
    x1 = arange(a+h,b,h)
    return zip(x0,x1)
    

#
# FEA2D
#
# SYNOPSIS
# One dimensional finite element analysis of the second order constant-coefficient
# equation with Dirichlet boundary conditions.
#
#
# INPUT
# el    - Intervals (pairs) partitioning some interval.
# p     - Element degrees.
# kappa - Constant coefficient.
# f     - Input function.
# xw    - Quadrature rules for evaluating right-hand side functionals.
#
# RETURNS
# Solution written in element bases. 
#
def fea_diri2(els, p, kappa, f, xw):

    # store element lengths
    hs = map(lambda el : el[1]-el[0], els)

    # reference element nodes
    xi = linspace(-1,1,p+1)
    # compute Lagrange shape functions and derivatives
    phi = lagrange(p+1)
    phid = matrix([polyder(phi[i,:].tolist()[0]) for i in range(0,p+1)])

    # determine local matrices
    loc0 = zeros(p+1)
    loc1 = zeros(p+1)
    for k,l in product(range(p+1),range(p+1)):
        loc0[k,l] = quadpq(phi[k,:], phi[l,:], xw)
        loc1[k,l] = quadpq(phid[k,:], phid[l,:], xw)
    
    print phi
    print phid
    print loc0
    print loc1

    dof = (len(els)+1)*(p+1)
    print "dof = ", dof, "nel =", len(els)

    nel = len(els)

    ind = map(lambda i : range(i*(p+1),(i+1)*(p+1)), range(nel))

    print ind

    A = zeros(dof)
    b = zeros(dof,1)

    # Assembly
    e = 0
    for el in els:
        Jaff = mpf(p)/hs[e]
        loc = Jaff*(kappa*loc0) + (1/Jaff)*loc1

        ia = ind[e][0]
        ib = ind[e][p]+1
        print ia,ib
        A[ia:ib,ia:ib] += loc

        # Deduce right hand functional
        for nu in range(p+1):
            ft = polyaff(f, -1, 1, el[0], el[1])
            b[e+nu] += Jaff*quadpq(ft,phi[nu,:],xw)

        # Dirichlet boundary conditions
        if e == 0:
            b[e+1] -= loc0[0,1]*polyval(f, el[0])
        if e+1 == len(els):
            b[e+0] -= loc0[1,0]*polyval(f, el[1])

        e += 1

    print A,b
    
    return
