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
