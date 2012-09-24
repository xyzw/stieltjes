from matrix_utils import *
from poly import *
from ppoly import *
from quad import *
from itertools import product
from fea1d import *
from bub import *
import pylab as pl
import copy
import time


def ppolympmathquad(pp,rhs,r0,r1,el0,el1):
    s = 0
    for i in range(len(pp.intv)):
        s += rhs(pp.poly[i],pp.intv[i][0],pp.intv[i][1],el0,el1)
    return s

def spikeder(s):
    spikeder = lambda x : s*(1-x**2)**(s-1)*(-2*x)
    return spikeder


# Solve implicit a posteriori Neumann problem
def iapneu(els,G,x,phi,kappa2,rhs,duX,iapnel,iapp):
    nel = len(els)
    uh = ppolyfea1sol(els,G,x,phi)
    duh = ppolyder(uh)
    res = ppolyaxpy(-kappa2,uh,ppolyder(duh)) # u_h''-kappa2 u_h

    xwl = gauss(r_jacobi(phi.cols+iapp))

    h1norm2 = []
    eh = ppoly()

    k = 0
    for el in els:
        psi = lagrangecheby(iapp)
 
        # Approximate derivatives on the boundaries
        if duX is None:
            A = 0 if k == 0 else 0.5*(polyval(duh.poly[k],el[0])-polyval(duh.poly[k-1],el[0]))
            B = 0 if k == nel-1 else 0.5*(polyval(duh.poly[k],el[1])-polyval(duh.poly[k+1],el[1]))
        else:
            A = duX[k]-polyval(duh.poly[k],el[0])
            B = duX[k+1]-polyval(duh.poly[k],el[1])

        #print "Ae={0:s}, Be={0:s}".format(nstr(Ae),nstr(Be))
        
        erhs = lambda phii, r0, r1, el0, el1 : rhs(phii, r0, r1, el0, el1) +\
                0.5*(el1-el0) * quadpq(phii,polyaff(res.poly[k],mpf(-1),mpf(1),el0,el1), xwl)

        Xs = linspace(el[0],el[1],iapnel+1)
        subels = zip(Xs[0:iapnel], Xs[1:iapnel+1])

        eels,eG,ex,epsi = fea1dapel1(subels, psi, kappa2, 1, [A, B], erhs)

        eeh = ppolyfea1sol(eels,eG,ex,epsi)
        h1norm2.append(ppolyh1norm2(eeh,xwl))
        
        ppolyext(eh, eeh)

        k += 1

    return eh, h1norm2

def iapbub(els,G,x,phi,kappa2,rhs,iapnel,iapp,sobolev=True):
    nel = len(els)
    uh = ppolyfea1sol(els,G,x,phi)
    duh = ppolyder(uh)
    res = ppolyaxpy(-kappa2,uh,ppolyder(duh)) # u_h''-kappa2 u_h
    p = phi.cols

    #print ">>>> iapbub on", els

    # Precompute bubble functions
    mp.dps = 20
    alpha = p/2+1

    h = (els[0][1]-els[0][0]) # ASSUME UNIFORM ELEMENTS
    nbub = iapp
    
    if sobolev == True:
        bub = bubble(nbub,h,kappa2,alpha,p,iapnel)
    else:
        bub = bubblexn(nbub,h,kappa2,alpha,p,iapnel)
        
    degbub = len(bub[0].poly[0])

    # Determine local matrix
    xw = gauss(r_jacobi(degbub))
    dbub = [ppolyder(b) for b in bub]
    loc = zeros(nbub)
    for k,l in product(range(nbub),range(nbub)):
        loc[k,l] = kappa2*quadppqq(bub[k], bub[l], xw) + quadppqq(dbub[k], dbub[l], xw)

    # for the quadrature of the residual
    xwl = gauss(r_jacobi(max(degbub,p)))

    #nprint(loc)

    mp.dps = 15

    h1norm2 = []
    eh = ppoly()

    t1 = time.time()

    k = 0
    for el in els:
        # Translate bubbles to the element
        tbub = []
        for b in bub:
            bt = ppolytrans(b, 0.5*(el[0]+el[1]))
            tbub.append(bt)

        # Execute bubble kernel
        erhs = lambda phii, el0, el1 : rhs(phii, el0, el1, el0, el1) +\
               0.5*(el1-el0) * quadpqab(phii,res.poly[k],el0,el1, xwl)
        
        G,x = iapbubker(tbub, loc, erhs)
        eeh = ppolyiapbubsol(G,x,tbub)
        h1norm2.append(ppolyh1norm2(eeh,xwl))
        ppolyext(eh, eeh)

        k += 1
                  
    return eh,h1norm2,time.time()-t1

def iapbuba(els,G,x,phi,kappa2,rhs,iapp):
    nel = len(els)
    uh = ppolyfea1sol(els,G,x,phi)
    duh = ppolyder(uh)
    res = ppolyaxpy(-kappa2,uh,ppolyder(duh)) # u_h''-kappa2 u_h
    p = phi.cols


    # Precompute bubble functions
    mp.dps = 30
    alpha = p/2+1

    h = (els[0][1]-els[0][0]) # ASSUME UNIFORM ELEMENTS
    nbub = iapp


    xwl = gauss(r_jacobi(nbub+2*alpha))
    P = prebub(nbub,h,kappa2,alpha,xwl)
    bub = []
    for i in range(P.rows):
        bub.append(polytoppoly(P[i,:], [(-0.5*h, 0.5*h)]))
        
    degbub = len(bub[0].poly[0])
    # for the quadrature of the residual
    xwl = gauss(r_jacobi(max(degbub,p)))

    #nprint(loc)

    mp.dps = 15

    h1norm2 = []
    eh = ppoly()

    t1 = time.time()

    k = 0
    for el in els:
        # Translate bubbles to the element
        tbub = []
        for b in bub:
            bt = ppolytrans(b, 0.5*(el[0]+el[1]))
            tbub.append(bt)

        # Execute bubble kernel
        erhs = lambda phii, el0, el1 : rhs(phii, el0, el1, el0, el1) +\
               0.5*(el1-el0) * quadpqab(phii,res.poly[k],el0,el1, xwl)
        
        G,x = iapbubaker(tbub, erhs)
        eeh = ppolyiapbubsol(G,x,tbub)
        h1norm2.append(ppolyh1norm2(eeh,xwl))
        ppolyext(eh, eeh)

        k += 1
                  
    return eh,h1norm2,time.time()-t1

def iapbubaker(bub, rhs):    
    nbub = len(bub)
    G = [range(nbub)]       
    x = zeros(nbub,1)

    # Solve
    for i in range(nbub):
        bi = 0
        for k in range(len(bub[i].intv)):
            bi += rhs(bub[i].poly[k], bub[i].intv[k][0], bub[i].intv[k][1])
        x[i] = bi

    return G,x

def iapbubortho(els,G,x,phi,kappa2,rhs,iapnel,iapp):
    nel = len(els)
    uh = ppolyfea1sol(els,G,x,phi)
    duh = ppolyder(uh)
    res = ppolyaxpy(-kappa2,uh,ppolyder(duh)) # u_h''-kappa2 u_h
    p = phi.cols


    # Precompute bubble functions
    mp.dps = 30
    alpha = p/2+1

    h = (els[0][1]-els[0][0]) # ASSUME UNIFORM ELEMENTS
    nbub = iapp

    bub,bubinner = bubbleortho(nbub,h,kappa2,alpha,p,iapnel)
  
    degbub = len(bub[0].poly[0])
    
    # for the quadrature of the residual
    xwl = gauss(r_jacobi(max(degbub,p)))

    mp.dps = 15

    h1norm2 = []
    eh = ppoly()

    t1 = time.time()

    k = 0
    for el in els:
        # Translate bubbles to the element
        tbub = []
        for b in bub:
            bt = ppolytrans(b, 0.5*(el[0]+el[1]))
            tbub.append(bt)

        # Execute bubble kernel
        erhs = lambda phii, el0, el1 : rhs(phii, el0, el1, el0, el1) +\
               0.5*(el1-el0) * quadpqab(phii,res.poly[k],el0,el1, xwl)
        
        G,x = iapbuborthoker(tbub, bubinner, erhs)
        eeh = ppolyiapbubsol(G,x,tbub)
        h1norm2.append(ppolyh1norm2(eeh,xwl))
        ppolyext(eh, eeh)

        k += 1
                  
    return eh,h1norm2,time.time()-t1

def iapbuborthoker(bub, bubinner, rhs):    
    nbub = len(bub)
    G = [range(nbub)]       
    x = zeros(nbub,1)

    # Solve
    for i in range(nbub):
        bi = 0
        for k in range(len(bub[i].intv)):
            bi += rhs(bub[i].poly[k], bub[i].intv[k][0], bub[i].intv[k][1])
        x[i] = bi/bubinner[i]

    return G,x


def iapbubxnortho(els,G,x,phi,kappa2,rhs,iapnel,iapp):
    nel = len(els)
    uh = ppolyfea1sol(els,G,x,phi)
    duh = ppolyder(uh)
    res = ppolyaxpy(-kappa2,uh,ppolyder(duh)) # u_h''-kappa2 u_h
    p = phi.cols


    # Precompute bubble functions
    mp.dps = 30
    alpha = p/2+1

    h = (els[0][1]-els[0][0]) # ASSUME UNIFORM ELEMENTS
    nbub = iapp

    bub,bubinner = bubblexnortho(nbub,h,kappa2,alpha,p,iapnel)
  
    degbub = len(bub[0].poly[0])
    
    # for the quadrature of the residual
    xwl = gauss(r_jacobi(max(degbub,p)))

    mp.dps = 15

    h1norm2 = []
    eh = ppoly()

    t1 = time.time()

    k = 0
    for el in els:
        # Translate bubbles to the element
        tbub = []
        for b in bub:
            bt = ppolytrans(b, 0.5*(el[0]+el[1]))
            tbub.append(bt)

        # Execute bubble kernel
        erhs = lambda phii, el0, el1 : rhs(phii, el0, el1, el0, el1) +\
               0.5*(el1-el0) * quadpqab(phii,res.poly[k],el0,el1, xwl)
        
        G,x = iapbuborthoker(tbub, bubinner, erhs)
        eeh = ppolyiapbubsol(G,x,tbub)
        h1norm2.append(ppolyh1norm2(eeh,xwl))
        ppolyext(eh, eeh)

        k += 1
                  
    return eh,h1norm2,time.time()-t1

def iapbubxnorthoker(bub, bubinner, rhs):    
    nbub = len(bub)
    G = [range(nbub)]       
    x = zeros(nbub,1)

    # Solve
    for i in range(nbub):
        bi = 0
        for k in range(len(bub[i].intv)):
            bi += rhs(bub[i].poly[k], bub[i].intv[k][0], bub[i].intv[k][1])
        x[i] = bi/bubinner[i]

    return G,x

def iapbubspread(els,G,x,phi,kappa2,rhs,iapnel,iapp,sobolev=True):
    nel = len(els)
    uh = ppolyfea1sol(els,G,x,phi)
    duh = ppolyder(uh)
    res = ppolyaxpy(-kappa2,uh,ppolyder(duh)) # u_h''-kappa2 u_h
    p = phi.cols

    #print ">>>> iapbub on", els

    # Precompute bubble functions
    mp.dps = 20
    alpha = p/2+1

    # ASSUME UNIFORM ELEMENTS !!!

    h1 = (els[0][1]-els[0][0]) # element size
    h = iapnel*h1  # spread size
    nbub = iapp
    spread = (iapnel-1)/2
    
    if sobolev == True:
        bub1 = bubble(nbub,h1,kappa2,alpha,p,iapnel)
        bub = bubble(nbub,h,kappa2,alpha,p,iapnel)
    else:
        bub1 = bubblexn(nbub,h1,kappa2,alpha,p,iapnel)
        bub = bubblexn(nbub,h,kappa2,alpha,p,iapnel)
        
    degbub = len(bub[0].poly[0])

    # Determine local matrices
    xw = gauss(r_jacobi(degbub))
    dbub1 = [ppolyder(b) for b in bub1]
    dbub = [ppolyder(b) for b in bub]
    loc1 = zeros(nbub)
    loc = zeros(nbub)
    for k,l in product(range(nbub),range(nbub)):
        loc1[k,l] = kappa2*quadppqq(bub1[k], bub1[l], xw) + quadppqq(dbub1[k], dbub1[l], xw)
        loc[k,l] = kappa2*quadppqq(bub[k], bub[l], xw) + quadppqq(dbub[k], dbub[l], xw)

    # for the quadrature of the residual
    xwl = gauss(r_jacobi(max(degbub,p)))

    #nprint(loc)

    mp.dps = 15

    h1norm2s = []
    eh = ppoly()

    t1 = time.time()

    k = 0
    for el in els:
        spreadels = []

        # Only use spreaded bubbles if the whole spread interval is in the domain
        
        if k-spread < 0 or k+spread >= nel:
            # Fall back to 0 spread case
            
            # Translate bubbles to the element
            tbub = []
            for b in bub1:
                bt = ppolytrans(b, 0.5*(el[0]+el[1]))
                tbub.append(bt)

            erhs = lambda phii, el0, el1 : rhs(phii, el0, el1, el0, el1) +\
                   0.5*(el1-el0) * quadpqab(phii,res.poly[k],el0,el1, xwl)
        
            G,x = iapbubker(tbub, loc1, erhs)
            eeh = ppolyiapbubsol(G,x,tbub)
            h1norm2s.append(ppolyh1norm2(eeh,xwl))
            ppolyext(eh, eeh)
        else:
            # Multiple element support
            tbub = []
            for b in bub:
                bt = ppolytrans(b, 0.5*(el[0]+el[1]))
                tbub.append(bt)
#                    xx,yy = ppolyvalres(bt, 100)
#                    pl.plot(xx,yy,label="spread",linewidth=1.5)

            erhs = lambda phii, l, el0, el1 : rhs(phii, el0, el1, el0, el1) +\
                   0.5*(el1-el0) * quadpqab(phii,res.poly[k-spread+l],el0,el1, xwl)
        
            G,x = iapbubspreadker(tbub, loc, erhs)
            eeh = ppolyiapbubsol(G,x,tbub)

            # We are only concerened with the restriction to the current element
            h1norm2s.append(h1norm2ab(eeh.poly[spread],eeh.intv[spread][0],eeh.intv[spread][1],xwl))
            eh.intv.append(eeh.intv[spread])
            eh.poly.append(eeh.poly[spread])
            #xx,yy = polyvalres(eeh.poly[spread], eeh.intv[spread][0], eeh.intv[spread][1], 1000)
            #pl.plot(xx,yy,label="spread",linewidth=2)
            #xx,yy = ppolyvalres(eeh, 1000)
            #pl.plot(xx,yy,label="spread",linewidth=2)

        k += 1
                  
    return eh,h1norm2s,time.time()-t1

# translate ppoly by a
def ppolytrans(pp, a):
    qq = ppoly()
    for i in range(len(pp.intv)):
        qq.intv.append((pp.intv[i][0]+a, pp.intv[i][1]+a))
        qq.poly.append(polycomp(pp.poly[i], poly([1, -a])))
    return qq

#def ppolyiapbubsol(els,G,x,phi):
#    pp = ppoly(els)
#    for e in range(len(els)):
#        q = poly([0])
#        for p in G[e]:
#            q = polyaxpy(x[p],phi[p,:],q)            
#        pp.poly[e] = q
#
 #   return pp

def ppolyiapbubsol(G,x,bub):
    pp = ppoly(bub[0].intv)
    for i in G[0]:
        pp = ppolyaxpy(x[G[0][i]], bub[i], pp)
    return pp

def iapbubspreadker(bub, loc, rhs):
    #print ">>>> iapbubker on [{0:s},{1:s}]".format(x0,x1)
    
    nbub = len(bub)

    # Generate Local-to-Global index map
    G = [range(nbub)]       
    dof = nbub
        
    #print els
    #print G

    A = zeros(dof)
    b = zeros(dof,1)

    # Assembly
    for i in range(nbub):
        for j in range(nbub):
            A[G[0][i],G[0][j]] += loc[i,j]

        for k in range(len(bub[i].intv)):
            b[G[0][i]] += rhs(bub[i].poly[k], k, bub[i].intv[k][0], bub[i].intv[k][1])

    #print ">>>> A"
    #print nstr(chop(A))
    #print ">>>> b"
    #print nstr(chop(b))

    x = lu_solve(A,b)
    
    #print ">>>> x"
    #print x
    
    return G,x



def iapbubker(bub, loc, rhs):
    #print ">>>> iapbubker on [{0:s},{1:s}]".format(x0,x1)
    
    nbub = len(bub)

    # Generate Local-to-Global index map
    G = [range(nbub)]       
    dof = nbub
        
    #print els
    #print G

    A = zeros(dof)
    b = zeros(dof,1)

    # Assembly
    for i in range(nbub):
        for j in range(nbub):
            A[G[0][i],G[0][j]] += loc[i,j]

        for k in range(len(bub[i].intv)):
            b[G[0][i]] += rhs(bub[i].poly[k], bub[i].intv[k][0], bub[i].intv[k][1])

    #print ">>>> A"
    #print nstr(chop(A))
    #print ">>>> b"
    #print nstr(chop(b))

    x = lu_solve(A,b)
    
    #print ">>>> x"
    #print x
    
    return G,x

def fea1dapel2(els, psi, kappa2, bt, d, rhs):
    n = len(els)+1
    nel = len(els)
    npsi = psi.rows/nel # basis functions per element
    degpsi = psi.cols # basis degree
    hs = map(lambda el : el[1]-el[0], els)

    # Element endpoints
    x0 = els[0][0]
    x1 = els[-1][1]

    #print ">>>> [{0:s},{1:s}] -- executing bubble a posteriori kernel".format(nstr(x0),nstr(x1))
    #print nstr(psi)

    # Determine local matrices for every subelement's basis
    xw = gauss(r_jacobi(degpsi))

    # Differentiate basis
    dpsi = zeros(psi.rows,psi.cols-1)
    for i in range(psi.rows):
        dpsi[i,:] = polyder(psi[i,:])

    # Compute local matrices for every element
    loc0 = [zeros(npsi) for i in range(nel)]
    loc1 = [zeros(npsi) for i in range(nel)]
    
    for i in range(nel):
        for k,l in product(range(npsi),range(npsi)):
            el0 = els[i][0]
            el1 = els[i][1]
            psik = polyaff(psi[i*npsi+k,:], mpf(-1), mpf(1), el0, el1)
            psil = polyaff(psi[i*npsi+l,:], mpf(-1), mpf(1), el0, el1)
            dpsik = polyaff(dpsi[i*npsi+k,:], mpf(-1), mpf(1), el0, el1)
            dpsil = polyaff(dpsi[i*npsi+l,:], mpf(-1), mpf(1), el0, el1)
            loc0[i][k,l] = quadpq(psik, psil, xw)
            loc1[i][k,l] = quadpq(dpsik, dpsil, xw)

    #print ">>>> phi"
    #print phi
    #print phid
    #print ">>>> loc0"
    print nstr(loc0)
    print ">>>> loc1"
    print nstr(loc1)

    # Generate Local-to-Global index map
    G = []

    if bt == 0:
        #G.extend([[i,i+1] for i in range(-1,nel-1)])
        #G[nel-1][1] = -1

        #for k in range(nel):
        #    G[k].extend(range(nel-1+k*(p-1),nel-1+(k+1)*(p-1)))

#        dof = nel-2 + nel*(p-1)

        for k in range(nel):
            G.append(range(k*npsi,(k+1)*npsi))
        
        dof = nel*npsi
    elif bt == 1:
        G.extend([[i,i+1] for i in range(nel)])

        for k in range(nel):
            G[k].extend(range(nel+1+k*(p-1),nel+1+(k+1)*(p-1)))
                
        dof = nel + nel*(p-1)
        
    #print els
    print ">>>> G"
    print G

    A = zeros(dof)
    b = zeros(dof,1)

    # Assembly
    for e in range(nel):
        Jaff = mpf(2)/hs[e]
        #loce = (mpf(1)/Jaff)*(kappa2*loc0[e]) + Jaff*loc1[e]
        loce = kappa2*loc0[e] + loc1[e]
        el0,el1 = els[e][0],els[e][1]
        print ">>>> ", el0, el1
        
        for i in range(npsi):
            for j in range(npsi):
                A[G[e][i],G[e][j]] += loce[i,j]
                    
            b[G[e][i]] += rhs(psi[npsi*e+i,:],el0,el1,el0,el1)

        ## elif bt == 1: # Neumann
        ##     if e == 0:
        ##         for k in range(0,p+1):
        ##             b[G[e][k]] -= polyval(bub[k].poly[e],x0)*d[0]
        ##     if e == nel-1:
        ##         for k in range(0,p+1):
        ##             b[G[e][k]] += polyval(bub[k].poly[e],x1)*d[1]

    print ">>>> A"
    print nstr(A)
    print ">>>> b"
    print nstr(b)

    x = lu_solve(A,b)
    
    ## if bt == 0:
    ##     x = col_join(x, matrix(d))
    ##     G[0][0] = dof+1
    ##     G[nel-1][p] = dof+2
    
    print ">>>> x"
    print x
    
    return els,G,x,psi

def fea1dapel1(els, phi, kappa2, bt, d, rhs):
    n = len(els)+1
    p = phi.rows-1
    hs = map(lambda el : el[1]-el[0], els)

    phid = matrix([polyder(phi[i,:].tolist()[0]) for i in range(0,p+1)])

    # determine local matrices if neccessary
    xw = gauss(r_jacobi(p+1))
    loc0 = zeros(p+1)
    loc1 = zeros(p+1)
    for k,l in product(range(p+1),range(p+1)):
        loc0[k,l] = quadpq(phi[k,:], phi[l,:], xw)
        loc1[k,l] = quadpq(phid[k,:], phid[l,:], xw)

    #print ">>>> phi"
    #print phi
    #print phid
    #print ">>>> loc0"
    #print loc0
    #print ">>>> loc1"
    #print loc1

    nel = len(els)

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
                b[G[e][i]] += rhs(phi[i,:],-1.0,1.0,el[0],el[1])

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


def fea1dapel0(X, phi, kappa2, bt, d, rhs):
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
