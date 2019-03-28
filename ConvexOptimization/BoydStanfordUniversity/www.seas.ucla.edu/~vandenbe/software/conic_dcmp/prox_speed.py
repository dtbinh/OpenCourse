from cvxopt import solvers
from matplotlib.pyplot import *
from proxlib import proxqp_clique
import matfile

try:
    import psutil
except(ImportError):
    assert False, "Python package psutil required. Package can be downloaded at http://code.google.com/p/psutil/"



"""
Collects runtime of prox evaluations using fast prox method for dense cliques. 
Solves

    minimize trace(C*X) + 1/2 * ||X-Z||_F^2
    subject to
        A(X) == b
        X >= 0
        
        
Requires problems to be first generated (via generate_single_prox.py) to maintain 
consistency with SDPT3 results.
"""

vec = range(25,251,25)
times = [-1. for i in vec]
obj = [-1. for i in vec]

for pin in xrange(len(vec)):
    p = vec[pin]
    s = p
    
    d = matfile.read("problems/p%d_s%d.mat" % (p,s))
    
    try:
        C, A, b, Z = d['C'],d['A'],d['b'],d['Z']
    except:
        assert False, "Problem p = %d, s = %d not yet generated. Please run generate_single_prox.py." % (p,s)
        
    start = psutil.cpu_times()[0]
    sol,obj[pin] = proxqp_clique(C, A, b, Z, 1.)
    times[pin] = psutil.cpu_times()[0] - start


figure(1)
clf()
plot(vec,times,'bo-')
xlabel('Order of clique')
ylabel('Runtime')
show()
  
