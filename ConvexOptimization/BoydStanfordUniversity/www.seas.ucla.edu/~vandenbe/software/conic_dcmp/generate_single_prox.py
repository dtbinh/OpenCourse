from cvxopt import matrix, normal
import matfile
import os
"""
Generates matrices A,b,C, Z for solving

    minimize trace(C*X) + 1/2 * ||X-Z||_F^2
    subject to
        A(X) == b
        X >= 0
        
and stores in problem directory. Run this before running prox_speed.py or
prox_speed.m

"""
os.mkdir('problems')
vec = range(25,251,25)

for pin in xrange(len(vec)):

    p = vec[pin]
    s = p
    print "generating p = %d, s = %d" % (p,s)    
    A = normal(p**2,s)
    for k in xrange(s):
        tmp = matrix(A[:,k],(p,p))
        tmp = tmp + tmp.T
        A[:,k] = tmp[:]
    
    Z = matrix(normal(p**2),(p,p))
    Z = Z*Z.T

    S = matrix(normal(p**2),(p,p))
    S = S*S.T
    
    b = A.T*Z[:]
    C = matrix(A*normal(s) + S[:],(p,p))
    
    matfile.write("problems/p%d_s%d.mat" % (p,s),{'C':C,'A':A,'b':b,'Z':Z})
    

