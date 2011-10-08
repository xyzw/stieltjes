from matrix_utils import *

#
# TQL2
#
# SYNOPSIS
# Computes the eigendecomposition of a symmetric tridiagonal matrix using
# the QL algorithm.
#
# INPUT
# d, e - Diagonal and subdiagonal as nx1 and (n-1)x1 vectors
# ev - Generate eigenvectors?
#
# OUTPUT
# d - Eigenvalues
#
# RETURN VALUE
# Matrix whose columns are eigenvectors.
#
# COMMENTS
# This is derived from the Algol procedures tql2, by
# Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
# Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
# Fortran subroutine in EISPACK.
#
def tql2(d,e,ev=False,maxit=50):
    if d.rows != e.rows+1:
        raise ValueError, "Invalid input vector dimension"
    
    n = d.rows
    e = col_join(e,zeros(1))
    V = zeros(n)

    f = mpf(0)
    tst1 = mpf(0)
    EPS = power(mpf(2), -mp.prec)

    for l in range(n):
        # Find small subdiagonal element
        tst1 = max(tst1,fabs(d[l]) + fabs(e[l]))
        m = l
        while m < n:
            if fabs(e[m]) <= EPS*tst1:
                break
            m += 1
   
         # If m == l, d[l] is an eigenvalue,
         # otherwise, iterate.
        if m > l:
            it = 0;
            while True:
                it = it + 1
                if it >= maxit:
                    raise ValueError, "Maximum iteration count reached without convergence"
   
                # Compute implicit shift
                g = d[l]
                p = (d[l+1] - g) / (mpf(2) * e[l])
                r = hypot(p,mpf(1))
                if p < mpf(0):
                    r = -r

                d[l] = e[l] / (p + r)
                d[l+1] = e[l] * (p + r)
                dl1 = d[l+1]
                h = g - d[l]
                for i in range(l+2,n):
                    d[i] -= h
               
                f = f + h
   
                # Implicit QL transformation.
                p = d[m]
                c = c2 = c3 = mpf(1)
                el1 = e[l+1]
                s = s2 = mpf(0)
                for i in range(m-1, l-1, -1):
                    c3 = c2
                    c2 = c
                    s2 = s
                    g = c * e[i]
                    h = c * p
                    r = hypot(p,e[i])
                    e[i+1] = s * r
                    s = e[i] / r
                    c = p / r
                    p = c * d[i] - s * g
                    d[i+1] = h + s * (c * g + s * d[i])
   
                    # Accumulate transformation.
                    for k in range(n):
                        h = V[k,i+1]
                        V[k,i+1] = s * V[k,i] + c * h
                        V[k,i] = c * V[k,i] - s * h
                        
                p = -s * s2 * c3 * el1 * e[l] / dl1
                e[l] = s * p
                d[l] = c * p
                
                # Check for convergence.
                if not fabs(e[l]) > EPS*tst1: break
         
        d[l] = d[l] + f
        e[l] = mpf(0)
     
    # Sort eigenvalues and corresponding vectors.
    for i in range(n-1):
        k = i
        p = d[i]
        for j in range(i+1, n):
            if d[j] < p:
                k = j
                p = d[j]
               
        if k != i:
            d[k] = d[i];
            d[i] = p;
            for j in range(n):
                p = V[j,i];
                V[j,i] = V[j,k];
                V[j,k] = p;

    return V

#
# SVD
#
# SYNOPSIS
# Computes the singular value decomposition of a square matrix.
#
# INPUT
# a - A square matrix
#
# RETURN VALUE
# u - Left eigenvector matrix
# q - Matrix of singular values
# v - Right eigenvector matrix
#
# COMMENTS
#
# Almost exact translation of the ALGOL SVD algorithm published in
# Numer. Math. 14, 403-420 (1970) by G. H. Golub and C. Reinsch
#
# Copyright (c) 2005 by Thomas R. Metcalf, helicity314-stitch <at> yahoo <dot> com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
import copy

def svd(a):
    # Golub and Reinsch state that eps should not be smaller than the
    # machine precision, ie the smallest number
    # for which 1+e>1.  tol should be beta/e where beta is the smallest
    # positive number representable in the computer.
    eps = mp.eps
    tol = power(10,-mp.prec+1)/mp.eps
    assert mpf(1.0)+eps > mpf(1.0) # if this fails, make eps bigger
    assert tol > mpf(0.0)     # if this fails, make tol bigger
    itmax = 100
    u = copy.deepcopy(a)
    m = a.rows
    n = a.cols
    #if __debug__: print 'a is ',m,' by ',n

    if m < n:
        if __debug__: print 'Error: m is less than n'
        raise ValueError,'SVD Error: m is less than n.'

    e = zeros(1,n)  # allocate arrays
    q = zeros(1,n)
    v = zeros(n,n)
 
    # Householder's reduction to bidiagonal form

    g = mpf(0)
    x = mpf(0)

    for i in range(n):
        e[i] = g
        s = mpf(0)
        l = i+1
        for j in range(i,m):
            s += (u[j,i]*u[j,i])
            
        if s <= tol:
            g = mpf(0)
        else:
            f = u[i,i]
            
            g = sqrt(s) if f < mpf(0) else -sqrt(s)
                
            h = f*g-s
            u[i,i] = f-g
            for j in range(l,n):
                s = mpf(0)
                for k in range(i,m):
                    s += u[k,i]*u[k,j]
                f = s/h
                for k in range(i,m):
                    u[k,j] = u[k,j] + f*u[k,i]

        q[i] = g
        s = mpf(0)

        for j in range(l,n):
            s = s + u[i,j]*u[i,j]
            
        if s <= tol:
            g = mpf(0)
        else:
            f = u[i,i+1]
            g = sqrt(s) if f < mpf(0) else -sqrt(s)
            h = f*g - s
            u[i,i+1] = f-g
            for j in range(l,n):
                e[j] = u[i,j]/h
                
            for j in range(l,m):
                s=mpf(0)
                for k in range(l,n):
                    s = s + u[j,k]*u[i,k]
                for k in range(l,n):
                    u[j,k] = u[j,k] + s*e[k]
                    
        y = fabs(q[i])+fabs(e[i])
        if y > x: x=y
    # accumulation of right hand transformations
    for i in range(n-1,-1,-1):
        if g != mpf(0):
            h = g*u[i,i+1]
            for j in range(l,n): v[j,i] = u[i,j]/h
            for j in range(l,n):
                s = mpf(0)
                for k in range(l,n): s += u[i,k]*v[k,j]
                for k in range(l,n): v[k,j] += s*v[k,i]
        for j in range(l,n):
            v[i,j] = mpf(0)
            v[j,i] = mpf(0)
        v[i,i] = mpf(1)
        g = e[i]
        l = i
    #accumulation of left hand transformations
    for i in range(n-1,-1,-1):
        l = i+1
        g = q[i]
        for j in range(l,n): u[i,j] = mpf(0)
        if g != mpf(0):
            h = u[i,i]*g
            for j in range(l,n):
                s=mpf(0)
                for k in range(l,m): s += u[k,i]*u[k,j]
                f = s/h
                for k in range(i,m): u[k,j] += f*u[k,i]
            for j in range(i,m): u[j,i] = u[j,i]/g
        else:
            for j in range(i,m): u[j,i] = mpf(0.0)
        u[i,i] += mpf(1.0)
    #diagonalization of the bidiagonal form
    eps = eps*x
    for k in range(n-1,-1,-1):
        for iteration in range(itmax):
            # test f splitting
            for l in range(k,-1,-1):
                goto_test_f_convergence = False
                if fabs(e[l]) <= eps:
                    # goto test f convergence
                    goto_test_f_convergence = True
                    break  # break out of l loop
                if fabs(q[l-1]) <= eps:
                    # goto cancellation
                    break  # break out of l loop
            if not goto_test_f_convergence:
                #cancellation of e[l] if l>0
                c = mpf(0)
                s = mpf(1)
                l1 = l-1
                for i in range(l,k+1):
                    f = s*e[i]
                    e[i] = c*e[i]
                    if abs(f) <= eps:
                        #goto test f convergence
                        break
                    g = q[i]
                    h = pythag(f,g)
                    q[i] = h
                    c = g/h
                    s = -f/h
                    for j in range(m):
                        y = u[j,l1]
                        z = u[j,i]
                        u[j,l1] = y*c+z*s
                        u[j,i] = -y*s+z*c
            # test f convergence
            z = q[k]
            if l == k:
                # convergence
                if z<mpf(0):
                    #q[k] is made non-negative
                    q[k] = -z
                    for j in range(n):
                        v[j,k] = -v[j,k]
                break  # break out of iteration loop and move on to next k value
            if iteration >= itmax-1:
                if __debug__: print 'Error: no convergence.'
                # should this move on the the next k or exit with error??
                #raise ValueError,'SVD Error: No convergence.'  # exit the program with error
                break  # break out of iteration loop and move on to next k
            # shift from bottom 2x2 minor
            x = q[l]
            y = q[k-1]
            g = e[k-1]
            h = e[k]
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(mpf(2.0)*h*y)
            g = pythag(f,mpf(1))
            if f < 0:
                f = ((x-z)*(x+z)+h*(y/(f-g)-h))/x
            else:
                f = ((x-z)*(x+z)+h*(y/(f+g)-h))/x
            # next QR transformation
            c = mpf(1)
            s = mpf(1)
            for i in range(l+1,k+1):
                g = e[i]
                y = q[i]
                h = s*g
                g = c*g
                z = pythag(f,h)
                e[i-1] = z
                c = f/z
                s = h/z
                f = x*c+g*s
                g = -x*s+g*c
                h = y*s
                y = y*c
                for j in range(n):
                    x = v[j,i-1]
                    z = v[j,i]
                    v[j,i-1] = x*c+z*s
                    v[j,i] = -x*s+z*c
                z = pythag(f,h)
                q[i-1] = z
                c = f/z
                s = h/z
                f = c*g+s*y
                x = -s*g+c*y
                for j in range(m):
                    y = u[j,i-1]
                    z = u[j,i]
                    u[j,i-1] = y*c+z*s
                    u[j,i] = -y*s+z*c
            e[l] = mpf(0)
            e[k] = f
            q[k] = x
            # goto test f splitting
        
            
    #vt = transpose(v)
    #return (u,q,vt)
    return (u,q,v)
