from matrix_utils import *
from poly import *
from ppoly import *
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
# X     - Nodes.
# p     - Element degrees.
# kappa - Constant coefficient.
# f     - Input function.
# xw    - Quadrature rules for evaluating right-hand side functionals.
#
# RETURNS
# Solution written in element bases. 
#
def fea_diri2(X, p, kappa, f, xw):
    n = len(X)
    els = zip(X[0:n-1], X[1:n])

    # store element lengths
    hs = map(lambda el : el[1]-el[0], els)

    # reference element nodes
    xi = linspace(-1,1,p+1)
    # compute Lagrange shape functions and derivatives
    Xc = chebyx(p+1)
    phi = lagrange(Xc)
    phid = matrix([polyder(phi[i,:].tolist()[0]) for i in range(0,p+1)])

    # determine local matrices
    loc0 = zeros(p+1)
    loc1 = zeros(p+1)
    for k,l in product(range(p+1),range(p+1)):
        loc0[k,l] = quadpq(phi[k,:], phi[l,:], xw)
        loc1[k,l] = quadpq(phid[k,:], phid[l,:], xw)
    
    #print phi
    #print phid
    #print loc0
    #print loc1

    nel = len(X)-1

    # Generate Local-to-Global index map
    G = []
    i = 0
    io = nel-1

    G = [[i] for i in range(io)]
    G.append([-1])
    for k in range(nel):
        G[k].extend(range(io+k*(p-1),io+(k+1)*(p-1)))
        G[k].append(k-1)
    
    dof = nel-2 + nel*(p-1)

    #print els
    print G

    A = zeros(dof+1)
    b = zeros(dof+1,1)

    # Assembly
    e = 0
    for el in els:
        Jaff = mpf(2*p)/hs[e]
        loc = Jaff*(kappa*loc0) + (1/Jaff)*loc1

        for i in range(p+1):
            for j in range(p+1):
                if G[e][i] != -1 and G[e][j] != -1:
                    A[G[e][i],G[e][j]] += loc[i,j]
                else:
                    if j == 0:
                        b[G[e][i]] -= (mpf(1)/Jaff)*loc0[0,i]*polyval(f, el[0])
                    else:
                        b[G[e][i]] -= (mpf(1)/Jaff)*loc0[p,i]*polyval(f, el[1])

            if G[e][i] != -1:
                ft = polyaff(f, -1, 1, el[0], el[1])
                b[G[e][i]] += Jaff*quadpq(ft,phi[i,:],xw)
        e += 1

    #print ">>>> A"
    #print A
   # print ">>>> b"
    #print b

    x = lu_solve(A,b)
    
   # print ">>>> x"
   # print x
    
    return els,G,x,phi

# Evaluate finite element solution x (DEPRECATED)
def evalfea1sol(els,G,x,phi,elres):
    xx = []
    yy = []
    e = 0
    print phi
    for el in els:
        xxx = linspace(el[0],el[1],elres)
        yyy = zeros(1,elres)
        for p in range(phi.rows):
            i = G[e][p]
            if i != -1:
                yyy += x[i]*polyvalv(polyaff(phi[p,:],el[0],el[1],-1,1),xxx)

        xx.extend(xxx)
        yy.extend(yyy)
        e += 1
    return xx,yy

# Convert fea1 solution to ppoly
def ppolyfea1sol(els,G,x,phi):
    pp = ppoly(els)
    e = 0
    for el in els:
        q = poly([0])
        for p in range(phi.rows):
            i = G[e][p]
            if i != -1:
                q = polyaxpy(x[i],polyaff(phi[p,:],el[0],el[1],-1,1),q)
        pp.poly[e] = q
        e += 1
    return pp
        
        
