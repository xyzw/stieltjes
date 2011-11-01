from matrix_utils import *

def chebyshev(N, mom, abm=None):
    if N<=0:
        raise ValueError, "N out of range"
    
    if N > mom/2:
        N = mom.cols/2
        
    if abm is None:
        abm = zeros(2*N-1,2)

    if N > (abm.rows+1)/2:
        N = (abm.rows+1)/2
 
    ab = zeros(N,2)
    ab[0,0] = abm[0,0]+mom[1]/mom[0]
    ab[0,1] = mom[0]

    normsq = zeros(N,1)
    if N == 1:
        normsq[0] = mom[0]
        return ab, normsq

    sig = zeros(N+1,2*N)
    sig[1,:] = mom[0:2*N]
    for n in range(3,N+2):
        for m in range(n-1,2*N-n+3):
            sig[n-1,m-1] = sig[n-2,m]-(ab[n-3,0]-abm[m-1,0])*sig[n-2,m-1]-ab[n-3,1]*sig[n-3,m-1]+abm[m-1,1]*sig[n-2,m-2]
        ab[n-2,0] = abm[n-2,0]+sig[n-1,n-1]/sig[n-1,n-2]-sig[n-2,n-2]/sig[n-2,n-3]
        ab[n-2,1] = sig[n-1,n-2]/sig[n-2,n-3]

    for n in range(1,N+1):
        normsq[n-1] = sig[n,n-1]

    return ab, normsq.T
