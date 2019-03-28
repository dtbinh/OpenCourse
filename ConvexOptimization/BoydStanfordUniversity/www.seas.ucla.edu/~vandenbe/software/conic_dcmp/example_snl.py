from cvxopt import matrix,solvers,spdiag, sparse, spmatrix, normal
from matplotlib.pyplot import *
from snl_problem import snl_problem ,get_unconverted_problem
from spingarn import solve 
from convert import convert
from misc import spdiag_rect
from proxlib import proxqp_clique_SNL
    
#############################################################################################################
#  Solve the Sensor Network Problem  for 100 sensors in 2 dimensions, with 3 neighbors and .1 noise level.  #
#############################################################################################################


#Problem parameters   
print "Generate problem"
n = 25

X, EDM, Spun = snl_problem(nsensors = n, noise=.2,dim = 3, neighbors = 3)

figure(1)
clf
hold(True)
for k in xrange(len(Spun.I)):
    i = Spun.I[k]
    j = Spun.J[k]
    plot(X[[i,j],0],X[[i,j],1],'k')
plot(X[:,0],X[:,1],'b*')   
title('Sensor network') 
print "(Close figure plots to continue)\n"
show()

#Generate problem
Aun,bun,Cun,dims_un = get_unconverted_problem(X, EDM, Spun)


#Solve unconverted problem using CVXOPT
print "Solve unconverted problem using CVXOPT"
usol = solvers.conelp(-matrix(bun), Aun, Cun[:], dims = dims_un)
print "Primal objective: %e" % -usol['primal objective']
print '(Press any key to continue)\n'
raw_input()

#Converted problem
print "Converted problem"
Acon,bcon,Ccon,Lcon, dims_con = convert(Aun,bun,Cun,Spun, dims_un, Atype = spmatrix)


#Solve convertd problem using CVXOPT
print "Solve converted problem using CVXOPT"
Acon2 = sparse([Acon.T,Lcon.T]).T
bcon2 = matrix([matrix(bcon),matrix(0.,(Lcon.size[1],1))])

csol = solvers.conelp(-bcon2, Acon2, Ccon, dims = dims_con)
print "Primal objective: %e" % -csol['primal objective']
print '(Press any key to continue)\n'
raw_input()

#Plot sparsity pattern for converted pattern.
print "Plot sparsity pattern for converted pattern."

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


#Paramters for spingarn's method
sigma = 1.0
rho = 1.75

#Solve problem using block diagonal correlative sparsity
print "Solve problem using block diagonal correlative sparsity"
sol1 = solve(Acon,bcon,Ccon,Lcon,dims_con,proxqp_clique_SNL,sigma,rho,reltol=1e-4,maxiter = 300, multiprocess = 4,blockdiagonal = True, log_cputime = False)
print "Primal objective: %e" % sol1['primal'][-1]
print "CPU time: %f, Walltime = %f" % (sol1['cputime'], sol1['walltime'])
print '(Press any key to continue)\n'
raw_input()

# Solve problem using block-diagonal correlative sparsity and adaptive step size
print "Solve problem using block arrow problem using block-diagonal correlative sparsity and adaptive step size"
sol2 = solve(Acon,bcon,Ccon,Lcon,dims_con,proxqp_clique_SNL,sigma,rho,reltol=1e-4,maxiter = 300, blockdiagonal = True, adaptive = True, multiprocess = 4, log_cputime = False)
print "Primal objective: %e" % sol2['primal'][-1]
print "CPU time: %f, Walltime = %f" % (sol2['cputime'], sol2['walltime'])
print '(Press any key to continue)\n'
raw_input()


# Plot results of Spingarn's method
print "Plot results of Spingarn's method"
figure(8)
clf()
subplot(3,1,1)
hold(True)
plot(sol1['primal'])
plot(sol2['primal'],'b--')
legend(['Steplength = 1.0','Steplength adaptive'])
xlabel('Iteration')
ylabel('Primal evolution')

subplot(3,1,2)
semilogy(sol1['rprimal'])
hold(True)
semilogy(sol1['rdual'],'r')
semilogy(sol2['rprimal'],'b--')
semilogy(sol2['rdual'],'r--')
xlabel('Iteration')
ylabel('Residual')
legend(['Primal (Steplength = 1.0)','Dual (Steplength = 1.0)', 'Primal (Steplength adaptive)','Dual (Steplength adaptive)'])

subplot(3,1,3)
plot(sol2['sigma'])
xlabel('Iteration')
ylabel('Step size')
print "(Close figure plots to continue)\n"
show()
