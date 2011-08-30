from mpmath import *
mp.dps = 50
print pi
print 2*asinh(j).imag
print gamma(0.5)**2
print sqrt(6*zeta(2))
print quad(lambda x: 4*sqrt(1-x**2), [0, 1])
print quad(lambda x: exp(-x**2), [-inf, inf]) ** 2
print nsum(lambda n: 4*(-1)**n/(2*n+1), [0, inf])
print limit(lambda n: 2**(4*n+1)*factorial(n)**4/(2*n+1)/factorial(2*n)**2, inf)
print findroot(sin, 3.14)
