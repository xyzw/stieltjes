from matrix_utils import *

# Compute recurrence coefficents for the Jacobi measure dmu(x)=(1-x)^a*(1+x)^b dx
def r_jacobi(N, a=0, b=0):    
    if (N <= 0) or (a <= -1) or (b <= -1):
        raise ValueError, "Exponents `a' and `b` must be strictly greater than negative one, and N must be nonnegative."
    
    nu = (b-a)/(a+b+2)
    mu = power(2,a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2)
    if N == 1:
        return [nu, mu]

    N = N-1
    n = matrix([linspace(1,N,N)])
    nu = matrix([nu])
    mu = matrix([mu])
    e = ones(1,N)
    nab = mpf(2)*n+(a+b)
    d = mult(nab,nab+2)
    A = row_join(nu, div((b*b-a*a)*e,d))

    n = n[1:N]
    nab = nab[1:N]
    
    B1=matrix([mpf(4)*(a+1)*(b+1)/((a+b+2)*(a+b+2)*(a+b+3))])
    B=mpf(4)*div(mult(mult(n+a,n+b),mult(n,n+a+b)), mult(mult(mult(nab,nab),nab+1), nab-1))
    ab=row_join(A.T, col_join(col_join(mu, B1),B.T));
    return ab
