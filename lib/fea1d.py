from matrix_utils import *
from poly import *
from ppoly import *
from quad import *
from itertools import product


# Returns a polynomial RHS lambda expression
rhspolyxw = None
def rhspoly(f,xw):
    return lambda phi,el0,el1 : ((el1-el0)/mpf(2)) * quadpq(polyaff(f, -1, 1, el0, el1),phi,xw)

def funaff(f, a, b, c, d):
    return lambda x : f(x*(c-d)/(a-b)+(a*d-b*c)/(a-b))

def rhsmpmathquad(f):
    return lambda phi,r0,r1,el0,el1 : (el1-el0)/(r1-r0) * quad(lambda x : funaff(f,r0,r1,el0,el1)(x)*polyval(phi, x), [r0, r1])

# Lagrange-Chebyshev basis polynomials of degree p
def lagrangecheby(p):
    Xc = cheby2(p)
    for i in range(p,1,-1):
        Xc[i],Xc[i-1] = Xc[i-1],Xc[i]
    phi = lagrange(Xc)
    return phi

#
# FEA1DH
#
# SYNOPSIS
# One dimensional finite element analysis of the second order constant-coefficient Helmholtz
# equation with Dirichlet/Neumann boundary conditions.
#
#
# INPUT
# X     - Nodes.
# phi   - Element basis functions
# kappa - Constant coefficient.
# bt    - 0=Dirichlet, 1=Neumann
# d     - Boundary values/derivatives
# rhs   - Input function: takes phi_i, el0, el1 and assumed to compute the integral phi_i(x)f(x)
#
# RETURNS
#
def fea1dh(X, phi, kappa2, bt, d, rhs):
    n = len(X)
    p = len(phi)-1
    els = zip(X[0:n-1], X[1:n])
    # store element lengths
    hs = map(lambda el : el[1]-el[0], els)

    phid = matrix([polyder(phi[i,:].tolist()[0]) for i in range(0,p+1)])

    # determine local matrices if neccessary
    xw = gauss(r_jacobi(p+1))
    loc0 = zeros(p+1)
    loc1 = zeros(p+1)
    for k,l in product(range(p+1),range(p+1)):
        loc0[k,l] = quadpq(phi[k,:], phi[l,:], xw)
        loc1[k,l] = quadpq(phid[k,:], phid[l,:], xw)

    #print d
    #print ">>>> phi"
    #print phi
    #print phid
    #print ">>>> loc0"
    #print loc0
    #print ">>>> loc1"
    #print loc1

    nel = len(X)-1

    # Generate Local-to-Global index map
    G = []

    if bt == 0:
        G.extend([[i,i+1] for i in range(-1,nel-1)])
        G[nel-1][1] = -1

        for k in range(nel):
            G[k].extend(range(nel-1+k*(p-1),nel-1+(k+1)*(p-1)))
        
        dof = nel-2 + nel*(p-1)
    elif bt == 1:
        G.extend([[i,i+1] for i in range(nel)])

        for k in range(nel):
            G[k].extend(range(nel+1+k*(p-1),nel+1+(k+1)*(p-1)))
                
        dof = nel + nel*(p-1)
        
    #print els
    #print G

    A = zeros(dof+1)
    b = zeros(dof+1,1)

    # Assembly
    e = 0
    for el in els:
        Jaff = mpf(2)/(el[1]-el[0])
        loc = (mpf(1)/Jaff)*(kappa2*loc0) + Jaff*loc1
        
        for i in range(p+1):
            for j in range(p+1):
                if G[e][i] != -1 and G[e][j] != -1:
                    A[G[e][i],G[e][j]] += loc[i,j]
                    
            if G[e][i] != -1:
                b[G[e][i]] += rhs(phi[i,:],mpf(-1),mpf(1),el[0],el[1])

        # Boundary conditions
        if bt == 0: # Dirichlet
            if e == 0:
                for k in range(1,p+1):
                    b[G[e][k]] -= loc[0,k]*d[0]
            elif e == nel-1:
                for k in [0] + range(2,p+1):
                    b[G[e][k]] -= loc[1,k]*d[1]
        elif bt == 1: # Neumann
            if e == 0:
                for k in range(0,p+1):
                    b[G[e][k]] -= polyval(phi[k,:],mpf(-1))*d[0]
            if e == nel-1:
                for k in range(0,p+1):
                    b[G[e][k]] += polyval(phi[k,:],mpf(1))*d[1]
                
        e += 1

    #print ">>>> A"
    #print chop(A)
    #print ">>>> b"
    #print b

    x = lu_solve(A,b)
    
    if bt == 0:
        x = col_join(x, matrix(d))
        G[0][0] = dof+1
        G[nel-1][1] = dof+2
    
   # print ">>>> x"
   # print x
    
    return els,G,x,phi

# Convert fea1 solution to ppoly
def ppolyfea1sol(els,G,x,phi):
    pp = ppoly(els)
    e = 0
    for el in els:
        q = zeros(1)
        for p in range(phi.rows):
            i = G[e][p]
            if i != -1:
                q = polyaxpy(x[i],polyaff(phi[p,:],el[0],el[1],-1,1),q)
        pp.poly[e] = q
        e += 1
    return pp
        
        
# Convert fea1 solution to ppoly
def ppolyfea1solpp(els,G,x,phi):
    pp = ppoly()
    r0 = phi[0].intv[0][0]
    r1 = phi[0].intv[-1][1]

    e = 0
    for el in els:
        qq = ppolyaff(ppolyzero(phi[0].intv), r0, r1, el[0], el[1])

        for p in range(len(phi)):
            if G[e][p] != -1:
                qq = ppolyaxpy(x[G[e][p]], ppolyaff(phi[p],el[0],el[1],r0,r1), qq)

        ppolyext(pp,qq)
        e += 1
    return pp
        
        
