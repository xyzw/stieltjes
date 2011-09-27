from matrix_utils import *

def bisection(a, b, epsilon, m1, m2):
    n = a.rows
    z = 0
    
    # extend the vector b with zeros so that b_0=0, b_n=0
    b = col_join(col_join(zeros(1), b), zeros(1))
    b2 = mult(b,b)
    
    # determine initial search interval I by Gershgorin's theorem
    f1 = lambda i : a[i]-(fabs(b[i])+fabs(b[i+1]))
    f2 = lambda i : a[i]+(fabs(b[i])+fabs(b[i+1]))
    
    xmin = min(f1(i) for i in range(n))
    xmax = max(f2(i) for i in range(n))

    # initial eigenvalue lower and upper bounds are xmin and xmax
    xlo = zeros(1,n)
    xhi = zeros(1,n)
    
    xhigh = xmax

    for i in range(m1,m2+1):
        xlo[i] = xmin
        xhi[i] = xmax

    # iterate on the eigenvalues backwards
    for k in range(m2,m1-1,-1):
        
        # determine a lower bound of the kth eigenvalue based on the previous bounds
        xlow = xmin
        for i in range(k,m1-1,-1):
            if xlo[i] > xlow:
                xlow = xlo[i]
                break

        if xhigh > xhi[k]:
            xhigh = xhi[k]

        # iterate while convergence
        while xhigh - xlow > 2*eps*(fabs(xlow) + fabs(xhigh)) + epsilon:
            
            # bisect
            xmid = (xlow+xhigh)/2

            # number of bisections
            z = z + 1

            # number of negative elements in the Sturm sequence
            c = -1

            # initial value
            q = 1
            
            # evaluate the Sturm sequence; the number c of negative elements
            # equals the number of eigenvalues smaller than xmid
            for i in range(n):
                # guard against division by zero
                q = a[i] - xmid - (b2[i]/q if not almosteq(q,0) else b2[i]/eps)
                if q < 0:
                    c = c + 1
                    
            if c < k:
                # xmid is an upper bound for the eigenvalues m1,m1+1,...,c and
                # an upper bound for eigenvalues c+1,c+2,...,k.
                xlow = xmid
                if c < m1:
                    xlo[m1] = xmid
                else:
                    xlo[c+1] = xmid
                    if xhi[c] > xmid:
                        xhi[c] = xmid
            else:
                # otherwise it is an upper bound
                xhigh = xmid

        xhi[k] = (xhigh + xlow)/2

    return xhi, z

