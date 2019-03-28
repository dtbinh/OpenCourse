from cvxopt import matrix, spmatrix, matrix, blas, amd, base, cholmod, misc, sparse, mul, solvers, sqrt
from proxlib import proxqp_blockdiagonal, proxqp_general
import time

         
def solve(A,b,C,L,dims,proxqp = None,sigma=1.0,rho=1.0, **kwargs):
    
    """
    
    Solves the SDP

        min.  < c, x > 
        s.t.  A(x) = b
              x >= 0

    and its dual

        max.  -< b, y > 
        s.t.  s >= 0.
             c + A'(y) = s 
    
    
    Input arguments.
    
        A   is an N x M sparse matrix where N = sum_i ns[i]**2 and M = sum_j ms[j]
            and ns and ms are the SDP variable sizes and constraint block lengths respectively.
            
            The expression A(x) = b can be written as A.T*xtilde = b, where
            xtilde is a stacked vector of vectorized versions of xi.
        
        b   is a stacked vector containing constraint vectors of 
                        size m_i x 1.
    
        C   is a stacked vector containing vectorized 'd' matrices 
            c_k of size n_k**2 x 1, representing symmetric matrices.



        L  is an N X P sparse matrix, where L.T*X = 0 represents the consistency
            constraints. If an index k appears in different cliques i,j, and
            in converted form are indexed by it, jt, then L[it,l] = 1, 
            L[jt,l] = -1 for some l.
            
        dims    is a dictionary containing conic dimensions.
            dims['l'] contains number of linear variables under nonnegativity constrant
            dims['q'] contains a list of quadratic cone orders (not implemented!)
            dims['s'] contains a list of semidefinite cone matrix orders
        
        proxqp   is either a function pointer to a prox implementation, or, if 
                the problem has block-diagonal correlative sparsity, a pointer 
                to the prox implementation of a single clique. The choices are:
                
                proxqp_general : solves prox for general sparsity pattern
                
                proxqp_clique : solves prox for a single dense clique with 
                                only semidefinite variables.
                
                proxqp_clique_SNL : solves prox for sensor network localization 
                                    problem
        
        sigma is a nonnegative constant (step size)
        
        rho is a nonnegative constaint between 0 and 2 (overrelaxation parameter)
        
        In addition, the following paramters are optional:
        
            maxiter : maximum number of iterations (default 100)
            
            reltol : relative tolerance (default 0.01). 
                        If rp < reltol and rd < reltol and iteration < maxiter, 
                        solver breaks and returns current value.
                        
            adaptive : boolean toggle on whether adaptive step size should be 
                        used. (default False)
            
            mu, tau, tauscale : parameters for adaptive step size (see paper)

            multiprocess : number of parallel processes (default 1). 
                            if multiprocess = 1, no parallelization is used.
                            
            blockdiagonal : boolean toggle on whether problem has block diagonal
                            correlative sparsity. Note that even if the problem
                            does have block-diagonal correlative sparsity, if
                            this parameter is set to False, then general mode 
                            is used. (default False)
                            

            verbose : toggle printout (default True)
            
            log_cputime : toggle whether cputime should be logged.
	    
	    
    The output is returned in a dictionary with the following files:
        
        x : primal variable in stacked form (X = [x0, ..., x_{N-1}]) where
            xk is the vectorized form of the nk x nk submatrix variable.
        
        y, z : iterates in Spingarn's method
        
        cputime, walltime : total cputime and walltime, respectively, spent in 
                            main loop. If log_cputime is False, then cputime is 
                            returned as 0.
        
        primal, rprimal, rdual : evolution of primal optimal value, primal 
                                residual, and dual residual (resp.)
        
        sigma : evolution of step size sigma (changes if adaptive step size is used.)
    

    """
    
    solvers.options['show_progress'] = False
    maxiter = kwargs.get('maxiter',100)
    reltol = kwargs.get('reltol',0.01)
    adaptive = kwargs.get('adaptive',False)
    mu = kwargs.get('mu',2.0)
    tau = kwargs.get('tau',1.5)
    multiprocess = kwargs.get('multiprocess',1)
    tauscale = kwargs.get('tauscale',0.9)
    blockdiagonal = kwargs.get('blockdiagonal',False)
    verbose = kwargs.get('verbose',True)
    log_cputime = kwargs.get('log_cputime',True)
    
    if log_cputime:
        try:
            import psutil
        except(ImportError): 
            assert False, "Python package psutil required to log cputime. Package can be downloaded at http://code.google.com/p/psutil/"
    
    #format variables
    nl,ns = dims['l'], dims['s']
    C = C[nl:]
    L = L[nl:,:]
    As, bs = [],[]
    cons = []
    offset = 0
    for k in xrange(len(ns)):
        Atmp = sparse(A[nl+offset:nl+offset+ns[k]**2,:])
        J = list(set(list(Atmp.J)))
        Atmp = Atmp[:,J]
        if len(sparse(Atmp).V) == Atmp[:].size[0]: Atmp = matrix(Atmp)
        else: Atmp = sparse(Atmp)
        As.append(Atmp)
        bs.append(b[J])
        cons.append(J)
        
        offset += ns[k]**2
    
    if blockdiagonal:
        if sum([len(c) for c in cons]) > len(b):
            print "Problem does not have block-diagonal correlative sparsity. Switching to general mode."
            blockdiagonal = False
            
    #If not block-diagonal correlative sprasity, represent A as a list of lists: 
    #   A[i][j] is a matrix (or spmatrix) if ith clique involves jth constraint block
    #Otherwise, A is a list of matrices, where A[i] involves the ith clique and 
    #ith constraint block only.
    
    if not blockdiagonal:
        while sum([len(c) for c in cons]) > len(b):
            tobreak = False
            for i in xrange(len(cons)):
                for j in xrange(i):
                    ci, cj = set(cons[i]), set(cons[j])
                    s1 = ci.intersection(cj)
                    if len(s1) > 0:
                        s2 = ci.difference(cj)
                        s3 = cj.difference(ci)
                        cons.append(list(s1))
                        if len(s2) > 0: 
                            s2 = list(s2)
                            if not (s2 in cons): cons.append(s2)
                        if len(s3) > 0:
                            s3 = list(s3)
                            if not (s3 in cons): cons.append(s3)
                        
                        cons.pop(i)
                        cons.pop(j)
                        tobreak = True

                        break
                if tobreak: break

       
        As,bs = [],[]
        for i in xrange(len(cons)):
            J = cons[i]
            bs.append(b[J])
            Acol = []
            offset = 0
            for k in xrange(len(ns)):
                Atmp = sparse(A[nl+offset:nl+offset+ns[k]**2,J])
                if len(Atmp.V) == 0: 
                    Acol.append(0)
                elif len(Atmp.V) == Atmp[:].size[0]: 
                    Acol.append(matrix(Atmp))
                else: 
                    Acol.append(Atmp)
                offset += ns[k]**2
            As.append(Acol)
        
    
    ms = [len(i) for i in bs]
    bs = matrix(bs)
    meq = L.size[1]
    
    if (not blockdiagonal) and multiprocess > 1: 
        print "Multiprocessing mode can only be used if correlative sparsity is block diagonal. Switching to sequential mode."
        multiprocess = 1
    
    assert rho > 0 and rho < 2, 'Overrelaxaton parameter (rho) must be (strictly) between 0 and 2'

    
    # create routine for projecting on { x | L*x = 0 }
    #{ x | L*x = 0 } -> P = I - L*(L.T*L)i *L.T
    LTL = spmatrix([],[],[],(meq,meq))
    offset = 0
    for k in ns:
        Lk = L[offset : offset + k**2,:]
        base.syrk(Lk, LTL, trans = 'T', beta=1.0)
        offset += k**2
    LTLi = cholmod.symbolic(LTL, amd.order(LTL))
    cholmod.numeric(LTL, LTLi)


    #y = y - L*LTLi*L.T*y
    nssq = sum(matrix([nsk**2 for nsk in ns]))
    def proj(y,ip=True):
        if not ip: y = +y
        tmp = matrix(0.0,size=(meq,1))   
        
        ypre = +y
        base.gemv(L,y,tmp,trans='T',\
            m = nssq, n = meq, beta = 1)
        
        cholmod.solve(LTLi, tmp)
        base.gemv(L,tmp,y,beta=1.0,alpha=-1.0,trans='N',\
            m = nssq, n = meq)
        if not ip: return y
    
    time_to_solve = 0
    
    #initialize variables
    X = C*0.0
    Y = +X
    Z = +X
    dualS = +X
    dualy = +b
    PXZ = +X
    
        
    proxargs = {'C':C,'A':As,'b':bs,'Z':Z,'X':X,'sigma':sigma,'dualS':dualS,'dualy':dualy,
            'ns':ns,'ms':ms,'multiprocess':multiprocess}
    
  
    if blockdiagonal: proxqp = proxqp_blockdiagonal(proxargs,proxqp)
    else: proxqp = proxqp_general

    if log_cputime: utime = psutil.cpu_times()[0]
    wtime = time.time()
    primal = []
    rpvec, rdvec = [],[]
    sigmavec = []
    for it in xrange(maxiter):
        pv, gap = proxqp(proxargs)
       
        blas.copy(Z,Y)
        blas.axpy(X,Y,alpha=-2.0)   
        proj(Y,ip = True)

        #PXZ = sigma*(X-Z)
        blas.copy(X,PXZ)
        blas.scal(sigma,PXZ)
        blas.axpy(Z,PXZ,alpha=-sigma)
        
        
        #z = z + rho*(y-x)
        blas.axpy(X,Y,alpha=1.0)
        blas.axpy(Y,Z,alpha=-rho)
        
       
        xzn = blas.nrm2(PXZ)
        xn = blas.nrm2(X)
        xyn = blas.nrm2(Y)
        proj(PXZ,ip = True)       

        rdual = blas.nrm2(PXZ)        
        rpri = sqrt(abs(xyn**2 - rdual**2))/sigma

        if log_cputime: cputime = psutil.cpu_times()[0] -utime
        else: cputime = 0
        
        walltime = time.time() -wtime
        
        if rpri / max(xn,1.0) < reltol and rdual / max(1.0,xzn) < reltol: break

      
        rpvec.append(rpri/max(xn,1.0))
        rdvec.append(rdual/max(1.0,xzn))
        primal.append(pv)
        if adaptive:
            if  (rdual/xzn*mu  < rpri/xn):
                sigmanew = sigma * tau
            elif (rpri/xn*mu < rdual/xzn):
                sigmanew = sigma / tau
            else: sigmanew = sigma
            if it % 10 == 0 and it > 0 and tau > 1.0:
                tauscale *= 0.9
                tau = 1 + (tau - 1)*tauscale
            sigma = max(min(sigmanew,10.0),0.1)
        sigmavec.append(sigma)
        if verbose: 
            if log_cputime: print "%d: primal = %e, gap = %e, (rp,rd) = (%e,%e), sigma = %f, (cputime,walltime) = (%f, %f)" % (it,pv,gap,rpri/max(xn,1.0),rdual/max(1.0,xzn),sigma,cputime, walltime)
            else: print "%d: primal = %e, gap = %e, (rp,rd) = (%e,%e), sigma = %f, walltime = %f" % (it,pv,gap,rpri/max(xn,1.0),rdual/max(1.0,xzn),sigma, walltime)

        
    sol = {}
    sol['x'] = X
    sol['y'] = Y
    sol['z'] = Z
    sol['cputime'] = cputime
    sol['walltime'] = walltime
    sol['primal'] = primal
    sol['rprimal'] = rpvec
    sol['rdual'] = rdvec
    sol['sigma'] = sigmavec
    return sol       
         
