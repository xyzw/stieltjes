from matrix_utils import *
from poly import *
from ppoly import *
from quad import *
from itertools import product
from fea1d import *
import pylab as pl

def ppolympmathquad(pp,rhs,r0,r1,el0,el1):
    s = 0
    for i in range(len(pp.intv)):
        s += rhs(pp.poly[i],pp.intv[i][0],pp.intv[i][1],el0,el1)
    return s

#
# FEA1DIAP
#
# 1d implicit a posteriori scheme 
#
# els - Elements
# G - Global index map
# x - solution vector
# phi - element basis
# bt - 0=Dirichlet, 1=Neumann
# d - Boundary values/derivatives
# rhs - RHS of the original problem presented as a lambda-expression of the form rhs(phi_i,el0,el1)
def fea1diap(els,G,x,iapnel,p,phi,kappa2,bt,d,rhs):
    nel = len(els)
    uh = ppolyfea1sol(els,G,x,phi)
    duh = ppolyder(uh)
    res = ppolyaxpy(-kappa2,uh,ppolyder(duh)) # u_h''-kappa2 u_h

    eh = ppoly()

    k = 0
    for el in els:
        print ">>>> Computing on ", el
        psi = []
        psia = lagrangecheby(p)
        for i in range(p+1):
            psi.append(polytoppoly(psia[i,:],[(-1,1)]))


        # Approximate derivatives on the boundaries
        if k != 0:
            A0 = polyval(duh.poly[k],el[0])
            A1 = polyval(duh.poly[k-1],el[0])
            print "duh|k(a)", A0
            print "duh|k-1(a)", A1
            A = 0.5*(A1-A0)
        else:
            A = -polyval(duh.poly[k],el[0])
            print "duh|k(a)", A

        if k+1 != nel:
            B0 = polyval(duh.poly[k],el[1])
            B1 = polyval(duh.poly[k+1],el[1])
            B = 0.5*(B1-B0)
            print "duh|k(b)", B0
            print "duh|k+1(b)", B1
        else:
            B = -polyval(duh.poly[k],el[1])
            print "duh|k(b)", B

        # Exhibit weak residual
        erhs = lambda phi, r0, r1, el0, el1 : ppolympmathquad(phi, rhs, r0, r1, el0, el1) +\
               ppolympmathquad(phi, rhsmpmathquad(lambda x : polyval(res.poly[k],x)), r0, r1, el0, el1) #+\
               #A*ppolyval(phi,-1.0) + B*ppolyval(phi,1.0)

        Xel = linspace(el[0],el[1],iapnel+1)
        eels,eG,ex,epsi = fea1dapel(Xel, psi, kappa2, 1, [A, B], erhs)
        
        eeh = ppolyfea1solpp(eels,eG,ex,epsi)
        ppolyext(eh, eeh)

        k += 1
                  
    return eh

def fea1dapel(X, phi, kappa2, bt, d, rhs):
    n = len(X)
    p = len(phi)-1
    els = zip(X[0:n-1], X[1:n])
    # store element lengths
    hs = map(lambda el : el[1]-el[0], els)

    phid = [ppolyder(phi[i]) for i in range(0,p+1)]

    # determine local matrices if neccessary
    xw = gauss(r_jacobi(p+1))
    loc0 = zeros(p+1)
    loc1 = zeros(p+1)
    for k,l in product(range(p+1),range(p+1)):
        loc0[k,l] = quadppqq(phi[k], phi[l], xw)
        loc1[k,l] = quadppqq(phid[k], phid[l], xw)

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
        Jaff = mpf(2)/hs[e]
        loc = (mpf(1)/Jaff)*(kappa2*loc0) + Jaff*loc1
        
        for i in range(p+1):
            for j in range(p+1):
                if G[e][i] != -1 and G[e][j] != -1:
                    A[G[e][i],G[e][j]] += loc[i,j]
                    
            if G[e][i] != -1:
                b[G[e][i]] += rhs(phi[i],-1.0,1.0,el[0],el[1])

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
                    b[G[e][k]] -= ppolyval(phi[k],mpf(-1))*d[0]
            if e == nel-1:
                for k in range(0,p+1):
                    b[G[e][k]] += ppolyval(phi[k],mpf(1))*d[1]

        e += 1

    #print ">>>> A"
    #print A
    #print ">>>> b"
    #print b

    x = lu_solve(A,b)
    
    if bt == 0:
        x = col_join(x, matrix(d))
        G[0][0] = dof+1
        G[nel-1][p] = dof+2
    
   # print ">>>> x"
   # print x
    
    return els,G,x,phi
