from cvxopt import matrix,solvers,spdiag, sparse, spmatrix, normal, setseed
from matplotlib.pyplot import *
from spingarn import solve 
from convert import convert
from misc import spdiag_rect
from proxlib import  proxqp_clique
from blockarrow_problem import blockarrow_problem

    
#######################################################################################################################
#  Solve a problem with block-arrow aggregate sparsity and block-diagonal correlative sparsity in the converted form. #
#######################################################################################################################


#Problem parameters   
p = 15   # size of diagonal blocks
r = 5   # no of cliques, vector varibles  (r+1 is no of block rows/columns)
s = 5   # no of constraints per clique
aw = 10   # arrow width ( = order of overlapped block)
print "Preparing to generate problem with %d diagonal blocks of size %d, %d constraints / clique, and %d arrow width." % (r,p,s,aw)

#Generate unconverted problem
print "Generate unconverted problem"
Aun,bun,Cun,Spun = blockarrow_problem(p,r,s,aw)
dims_un = {'l':0, 'q':[], 's': [p*r+aw]}
#Plot sparsity pattern of unconverted pattern. Schur complement system is dense
print "Plot sparsity pattern of unconverted pattern. Schur complement system is dense"
figure(1)
clf()
spy(matrix(Spun))
title('Aggregate sparsity')

figure(2)
clf()
spy(matrix(Cun))
title('Sparsity of C')

figure(3)
clf()
spy(matrix(Aun[:,10],Cun.size))
title('Sparsity of A_10')
print "(Close figure plots to continue)\n"
show()

#Solve unconverted problem using CVXOPT
print "Solve unconverted problem using CVXOPT"
solvers.options['show_progress'] = True
usol = solvers.conelp(-matrix(bun), Aun, matrix(Cun[:]), dims = dims_un)
print "Primal objective: %f" % -usol['primal objective']
print '(Press any key to continue)\n'
raw_input()

#Converted problem
print "Converted problem"
Acon,bcon,Ccon,Lcon, dims_con = convert(Aun,bun,Cun,Spun,dims_un,Atype = matrix)


#Plot sparsity pattern for converted pattern.
print "Plot sparsity pattern for converted pattern."

ctmp = []
offset = 0
for k in dims_con['s']:
    ctmp.append(matrix(Ccon[offset:offset + k**2],(k,k)))
    offset += k**2
figure(4)
clf()
spy(matrix(spdiag(ctmp)))
title('Sparsity of converted C')

figure(5)
clf()
spy(matrix(Acon))
title('Sparsity of constraints')

figure(6)
clf()
spy(matrix(Acon.T*Acon))
title('Correlative sparsity')
print "(Close figure plots to continue)\n"
show()


#Solve converted problem using CVXOPT
print "Solve converted problem using CVXOPT"
solvers.options['show_progress'] = True
csol = solvers.conelp(-matrix([bcon, matrix(0.,(Lcon.size[1],1))]), sparse([Acon.T,Lcon.T]).T, Ccon, dims = dims_con)
print "Primal objective: %e" % -csol['primal objective']
print '(Press any key to continue)\n'
raw_input()

#Paramters for spingarn's method
sigma = 1.0
rho = 1.75

# Solve problem using generic Spingarn solver
print "Solve problem using generic Spingarn solver"
sol1 = solve(Acon,bcon,Ccon,Lcon,dims_con,None,sigma,rho,reltol=1e-4,maxiter = 100, log_cputime = False)
print "Primal objective: %e" % sol1['primal'][-1]
print "CPU time: %f, Walltime = %f" % (sol1['cputime'], sol1['walltime'])
print '(Press any key to continue)'
raw_input()

#Solve problem using block diagonal correlative sparsity
print "Solve problem using block diagonal correlative sparsity"
sol2 = solve(Acon,bcon,Ccon,Lcon,dims_con,proxqp_clique,sigma,rho,reltol=1e-4,maxiter = 100, blockdiagonal = True, multiprocess = 4, log_cputime = False)
print "Primal objective: %e" % sol2['primal'][-1]
print "CPU time: %f, Walltime = %f" % (sol2['cputime'], sol2['walltime'])
print '(Press any key to continue)'
raw_input()

# Solve problem using block-diagonal correlative sparsity and adaptive step size
print "Solve problem using block arrow problem using block-diagonal correlative sparsity and adaptive step size"
sol3 = solve(Acon,bcon,Ccon,Lcon,dims_con,proxqp_clique,sigma,rho,reltol=1e-4,maxiter = 100, blockdiagonal = True, adaptive = True, multiprocess = 4, log_cputime = False)
print "Primal objective: %e" % sol3['primal'][-1]
print '(Press any key to continue)'
raw_input()


# Plot results of Spingarn's method
print "Plot results of Spingarn's method"
figure(8)
clf()
subplot(3,1,1)
hold(True)
plot(sol1['primal'])
plot(sol3['primal'],'b--')
legend(['Steplength = 1.0','Steplength adaptive'])
xlabel('Iteration')
ylabel('Primal evolution')

subplot(3,1,2)
semilogy(sol1['rprimal'])
hold(True)
semilogy(sol1['rdual'],'r')
semilogy(sol3['rprimal'],'b--')
semilogy(sol3['rdual'],'r--')
xlabel('Iteration')
ylabel('Residual')
legend(['Primal (Steplength = 1.0)','Dual (Steplength = 1.0)', 'Primal (Steplength adaptive)','Dual (Steplength adaptive)'])

subplot(3,1,3)
plot(sol3['sigma'])
xlabel('Iteration')
ylabel('Step size')
print "(Close figure plots to continue)"
show()
