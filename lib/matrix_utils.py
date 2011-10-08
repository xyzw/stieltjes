from mpmath import *

def row_join(A,B):
    if A.rows != B.rows:
        raise IndexError, "Input matrices with different row counts."
        
    C = matrix(A.rows, A.cols+B.cols)
    for i in range(A.rows): 
        for k in range(A.cols):
            C[i,k] = A[i,k]
        for k in range(B.cols):
            C[i,k+A.cols] = B[i,k]
    return C

def col_join(A,B):
    if A.cols != B.cols:
        raise IndexError, "Input matrices with different column counts."
    
    C = matrix(A.rows+B.rows, A.cols)
    for i in range(A.cols): 
        for k in range(A.rows):
            C[k,i] = A[k,i]
        for k in range(B.rows):
            C[k+A.rows,i] = B[k,i]
    return C

def mult(A,B):
    if A.cols != B.cols or A.rows != B.rows:
        raise IndexError, "Input matrices have different shapes."
    
    C = matrix(A.rows, A.cols)
    for i in range(A.rows): 
        for j in range(A.cols):
            C[i,j] = A[i,j]*B[i,j]
    return C

def div(A,B):
    if A.cols != B.cols or A.rows != B.rows:
        raise IndexError, "Input matrices have different shapes."
    
    C = matrix(A.rows, A.cols)
    for i in range(A.rows): 
        for j in range(A.cols):
            C[i,j] = A[i,j]/B[i,j]
    return C

def tridiag(a,b,c):
    if a.rows != b.rows+1 or b.rows != c.rows or a.rows != c.rows+1:
        raise ValueError, "Input vectors dimension mismatch"

    n = a.rows

    A = zeros(n,n)

    for i in range(n):
        A[i,i] = a[i]

    for i in range(n-1):
        A[i,i+1] = b[i]
        A[i+1,i] = c[i]

    return A

def jacobi_matrix(ab):
    if ab.cols != 2:
        raise ValueError, "Not a recurrence coefficient matrix specified as input"

    n = ab.rows
    sb = zeros(n-1,1)

    for i in range(1,n):
        sb[i-1] = sqrt(ab[i,1])

    return ab[:,0],sb
        
def hypot(a,b):
    if fabs(a) > fabs(b):
        r = b/a
        r = fabs(a)*sqrt(1+r*r)
    elif b != 0:
        r = a/b
        r = fabs(b)*sqrt(1+r*r)
    else:
        r = mpf(0)
    return r


def samesign(a,b):
    return fabs(a) if b >= mpf(0) else -fabs(a)
