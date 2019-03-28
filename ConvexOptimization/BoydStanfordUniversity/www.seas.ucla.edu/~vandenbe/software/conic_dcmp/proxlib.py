from cvxopt import matrix, blas, lapack, solvers, cholmod, spmatrix, printing, base,sparse,mul,div,sqrt, misc
from multiprocessing import Process,Queue
from misc import symmetrize, cngrnc, sgemv



def copy(x,y):
    #y := x
    if type(x) is matrix and type(y) is matrix:
        blas.copy(x,y)
    else:
        y[:] = 0.0
        y[x.I + x.size[0]*x.J] = x.V
def scal(alpha,x):
    #x := alpha*x
    if type(x) is matrix:
        blas.scal(alpha,x)
    else:
        x[x.I,x.J] = alpha*x.V

def proxqp_general(proxargs):
    """
    Solves the conic QP

        min.  < c, x > + (rho/2) || x - z ||^2
        s.t.  A(x) = b
              x >= 0

    and its dual

        max.  -< b, y > - 1/(2*rho) * || c + A'(y) - rho * z - s ||^2 
        s.t.  s >= 0.

    The primal variable is x = matrix(x_0[:], ...,  x_{N-1}[:]) with x_k a symmetric
    matrix of order nk.
    The dual variables are y = matrix(y_0, ..., y_{M-1}), with y_k a column
    vector of length mk, and s = matrix(s_0[:], .., s_{N-1}[:]), with s_i a symmetric
    matrix of order nk.

    In the primal cost function, c = (c_0[:], ..., c_{N-1}[:]), with c_k a 
    symmetric matrix, and < c, x > = sum_j tr (c_j * x_j).
    In the dual cost function, b = matrix(b_0, ..., b_{M-1}), with b_k a 
    column vector, and < b, y > = sum_j b_i' * y_i.

    The mapping A(x) is defined as follows. 
    The value of u = A(x) is u = matrix(u_0, ..., u_{M-1}) with u_i a vector
    defined by

        u_i = sum_j Aij(x_j),  i = 0, ..., M-1.

    The adjoint v = A'(y)  is v = (v_0, ..., v_{N-1}) with v_j a matrix
    defined as

        v_j = sum_i Aij'(y_i),  j = 0, ..., N-1.

    If we expand the primal and dual problems we therefore have

        min.  sum_j < c_j, x_j >  + (rho/2) * sum_j || x_j - z_j ||_F^2 
        s.t.  sum_j Aij(x_j) = b_i, i = 0, ..., M-1 
              x_j >= 0, j = 0, ..., N-1

    with variables x_j and its dual

        max.  sum_i b_i'*y_i - 1/(2*rho) *
              sum_j || c_j - sum_i Aij'(y_i) + rho *z_j - s_j ||_F^2
        s.t.  s_j >= 0, j = 0, ..., N-1

    with variables y_i, s_j.

    Input arguments.
    
        proxargs['C']   is a stacked vector containing vectorized 'd' matrices 
                        c_k of size n_k**2 x 1, representing symmetric matrices.

        proxargs['A']   is a list of a list of either 'd' matrices, 
            'd' sparse matrices, or 0.  
            
            If A[i][j] = 0, then variable block i 
            is not involved in constraint block j. 
            
            Otherwise, A[i][j] has size n_j**2 times m_i.  Each of its columns 
            represents a symmetric matrix of order n_j in unpacked column-major 
            order. The term  Aij(x[j]) in the primal constraint is given by

                Aij(x_j) = A[i][j]' * vec(x_j).

            The adjoint Aij'(y[i]) in the dual constraint is given by

                Aij'(y_i) = mat( A[i][j] * y_i ).


        proxargs['b']   is a stacked vector containing constraint vectors of 
                        size m_i x 1.

        proxargs['Z']   is a stacked vector containing variable Z vectors of 
                        size n_k**2 x 1.
            
        proxargs['sigma'] is a positive scalar (step size).
        
        proxargs['X']   contains the primal variable X in stacked vector form. 
        
        proxargs['dualy']   contains the dual variable y 
        
        proxargs['dualS']   contains the primal variable S
	
	    On output, proxargs['X'], proxargs['dualy'], and proxargs['dualZ'] 
	    will be filled with prox optimal primal and dual variables.
	    
    Output arguments.

        primal : trace(C*X) where X is the optimal variable
        
        gap : trace(X*S) as computed by  CVXOPT
	
    """
    
    c,A,b = proxargs['C'],proxargs['A'],proxargs['b']

    z,X = proxargs['Z'], proxargs['X']
    rho = proxargs['sigma']
    ns,ms = proxargs['ns'],proxargs['ms']
    multiprocess = proxargs['multiprocess']

    

    N = len(A[0])
    M = len(A)
    
    def Pf(u, v, alpha = 1.0, beta = 0.0):

        # v[k] := alpha * rho * u[k] + beta * v[k] 

        blas.scal(beta, v)
        blas.axpy(u, v, alpha = alpha*rho)
    q = (c - rho*z)[:]
   


   
    xp = +q[:]
    bz = +q[:]
    uy = +b[:]

    def Gf(u, v, alpha = 1.0, beta = 0.0, trans = 'N'):
        # v = -alpha*u + beta * v 
        blas.scal(beta, v)
        blas.axpy(u, v, alpha = -alpha)
        
    h = matrix(0.0, (sum(nk**2 for nk in ns), 1))
    dims = {'l': 0, 'q': [], 's': ns }
    
   
    def Af(u, v, alpha = 1.0, beta = 0.0, trans = 'N'):
        
        # v := alpha * A(u) + beta * v if trans is 'N'
        # v := alpha * A'(u) + beta * v if trans is 'T'
        blas.scal(beta, v)
        if trans == 'N': 
            offseti = 0
            for i in xrange(M):
                offsetj = 0
                for j in xrange(N):
                    if type(A[i][j]) is matrix or type(A[i][j]) is spmatrix:
                        sgemv(A[i][j], u, v, n = ns[j], m = ms[i],
                            trans = 'T', alpha = alpha, beta = 1.0,
                            offsetx = offsetj, offsety = offseti)
                    offsetj += ns[j]**2
                offseti += ms[i]
        else:
            offsetj = 0
            for j in xrange(N):
                offseti = 0
                for i in xrange(M):
                    if type(A[i][j]) is matrix or type(A[i][j]) is spmatrix:
                        sgemv(A[i][j], u, v, n = ns[j], m = ms[i],
                            trans = 'N', alpha = alpha, beta = 1.0,
                            offsetx = offseti, offsety = offsetj)
                    offseti += ms[i]
                offsetj += ns[j]**2

    def xdot(x, y): 

        offset = 0
        for k in xrange(N):
            misc.trisc(x, {'l': offset, 'q': [], 's': [ns[k]]})
            offset += ns[k]**2

        a = blas.dot(x,y)
        offset = 0
        for k in xrange(N):
            misc.triusc(x, {'l': offset, 'q': [], 's': [ns[k]]})
            symmetrize(x,n=ns[k],offset=offset)
            offset += ns[k]**2

        return a
    
    U = [ matrix(0.0, (nk,nk)) for nk in ns ]
    Vt = [ matrix(0.0, (nk,nk)) for nk in ns ]
    sv = [ matrix(0.0, (nk,1)) for nk in ns ]
    Gamma = [ matrix(0.0, (nk,nk)) for nk in ns ]
    As = [ [ +A[i][j] for j in xrange(N) ] for i in xrange(M) ]

    def F(W):

        for j in xrange(N):

            # SVD R[j] = U[j] * diag(sig[j]) * Vt[j]
            lapack.gesvd(+W['r'][j], sv[j], jobu = 'A', jobvt = 'A',
                U = U[j], Vt = Vt[j])               
            

            # Vt[j] := diag(sig[j])^-1 * Vt[j]
            for k in xrange(ns[j]):
                blas.tbsv(sv[j], Vt[j], n = ns[j], k = 0, ldA = 1,
                    offsetx = k*ns[j])
                    

            # Gamma[j] is an ns[j] x ns[j] symmetric matrix
            #
            #  (sig[j] * sig[j]') ./  sqrt(1 + rho * (sig[j] * sig[j]').^2)

            # S = sig[j] * sig[j]'
            S = matrix(0.0, (ns[j], ns[j]))
            blas.syrk(sv[j], S)
            Gamma[j][:] = div(S, sqrt(1.0 + rho * S**2))[:]
            symmetrize(Gamma[j], ns[j])


            # As represents the scaled mapping
            #
            #     As(x) = A(u * (Gamma .* x) * u')
            #    As'(y) = Gamma .* (u' * A'(y) * u)
            #
            # stored in a similar format as A, except that we use packed 
            # storage for the columns of As[i][j].
            
            for i in xrange(M):

                if (type(A[i][j]) is matrix )or ( type(A[i][j]) is spmatrix):
                
                    # As[i][j][:,k] = vec( 
                    #     (U[j]' * mat( A[i][j][:,k] ) * U[j]) .* Gamma[j])

                    copy(A[i][j], As[i][j])
                    As[i][j] = matrix(As[i][j])
                    for k in xrange(ms[i]):
                        cngrnc(U[j], As[i][j], trans = 'T', offsetx =
                            k * (ns[j]**2), n = ns[j])
                        blas.tbmv(Gamma[j], As[i][j], n = ns[j]**2,
                            k = 0, ldA = 1, offsetx = k * (ns[j]**2))

                    # pack As[i][j] in place
                    #pack_ip(As[i][j], ns[j])
                    for k in xrange(As[i][j].size[1]):
                        tmp = +As[i][j][:,k]
                        misc.pack2(tmp,{'l':0,'q':[],'s':[ns[j]]})
                        As[i][j][:,k] = tmp

                else:
                    As[i][j] = 0.0
                

        # H is an m times m block matrix with i, k block
        #
        #      Hik = sum_j As[i,j]' * As[k,j] 
        #
        # of size ms[i] x ms[k].  Hik = 0 if As[i,j] or As[k,j] 
        # are zero for all j

        H = spmatrix([], [], [], (sum(ms), sum(ms))) 
        rowid = 0  
        for i in xrange(M):
            colid = 0
            for k in xrange(i+1):
                sparse_block = True
                Hik = matrix(0.0, (ms[i], ms[k]))
                for j in xrange(N):
                    if (type(As[i][j]) is matrix) and \
                        (type(As[k][j]) is matrix):
                            sparse_block = False
                            # Hik += As[i,j]' * As[k,j]
                            if i == k:
                                blas.syrk(As[i][j], Hik, trans = 'T',
                                    beta = 1.0, k = ns[j]*(ns[j]+1)/2,
                                    ldA = ns[j]**2)
                            else:
                                blas.gemm(As[i][j], As[k][j], Hik,
                                    transA = 'T', beta = 1.0,
                                    k = ns[j]*(ns[j]+1)/2, 
                                    ldA = ns[j]**2, ldB = ns[j]**2)
                if not(sparse_block):
                    H[rowid:rowid+ms[i], colid:colid+ms[k]] \
                        = sparse(Hik)
                colid += ms[k]
            rowid += ms[i]    
        
        
        HF = cholmod.symbolic(H)
        cholmod.numeric(H, HF) 


        def solve(x, y, z): 
                      
            """
            Returns solution of 

                rho * ux + A'(uy) - r^-T * uz * r^-1 = bx
                A(ux)                                = by
                -ux               - r * uz * r'      = bz.

            On entry, x = bx, y = by, z = bz.
            On exit, x = ux, y = uy, z = uz.
            """


            
            # bz is a copy of z in the format of x
            blas.copy(z,bz)
            # x := x + rho * bz  
            #    = bx + rho * bz
            blas.axpy(bz, x, alpha = rho)

            # x := Gamma .* (u' * x * u)  
            #    = Gamma .* (u' * (bx + rho * bz) * u)
            offsetj = 0
            for j in xrange(N):
                cngrnc(U[j], x, trans = 'T', 
                    offsetx = offsetj, n = ns[j])
                blas.tbmv(Gamma[j], x, n = ns[j]**2, k = 0, 
                    ldA = 1, offsetx = offsetj) 
                offsetj += ns[j]**2
            
            # y := y - As(x)
            #   := by - As( Gamma .* u' * (bx + rho * bz) * u)
            
            blas.copy(x,xp)
            
             
            offsetj = 0
            for j in xrange(N):
                misc.pack2(xp,{'l':offsetj, 'q':[], 's':[ns[j]]})
                offsetj += ns[j]**2


            offseti = 0
            for i in xrange(M):
                offsetj = 0
                for j in xrange(N):
                    if type(As[i][j]) is matrix:
                        blas.gemv(As[i][j], xp, y, trans = 'T',
                            alpha = -1.0, beta = 1.0, 
                            m = ns[j]*(ns[j]+1)/2, n = ms[i],ldA = ns[j]**2, 
                            offsetx = offsetj, offsety = offseti) 
                    offsetj += ns[j]**2
                offseti += ms[i]
            # y := -y - A(bz) 
            #    = -by - A(bz) + As(Gamma .*  (u' * (bx + rho * bz) * u) 
            
            
            Af(bz, y, alpha = -1.0, beta = -1.0)  

            # y := H^-1 * y
            #    = H^-1 ( -by - A(bz) + As(Gamma.* u'*(bx + rho*bz)*u) )
            #    = uy 
            
            

            cholmod.solve(HF, y)


            
            # bz = Vt' * vz * Vt 
            #    = uz where
            # vz := Gamma .* ( As'(uy)  - x ) 
            #     = Gamma .* ( As'(uy)  - Gamma .* (u'*(bx + rho *bz)*u) )
            #     = Gamma.^2 .* ( u' * (A'(uy) - bx - rho * bz) * u ).
            blas.copy(x,xp)
            

            
            offsetj = 0
            for j in xrange(N):

                # xp is -x[j] = -Gamma .* (u' * (bx + rho*bz) * u)
                # in packed storage
                misc.pack2(xp,{'l':offsetj,'q':[],'s':[ns[j]]})
                offsetj += ns[j]**2    
            blas.scal(-1.0, xp)

            offsetj = 0
            for j in xrange(N):
                # xp +=  As'(uy)
                
                offseti = 0
                for i in xrange(M):
                    if type(As[i][j]) is matrix:
                        blas.gemv(As[i][j], y, xp, alpha = 1.0,
                             beta = 1.0, m = ns[j]*(ns[j]+1)/2, \
                                n = ms[i],ldA = ns[j]**2, \
                                offsetx = offseti, offsety = offsetj)
                    offseti += ms[i]

                
               
                # bz[j] is xp unpacked and multiplied with Gamma
                #unpack(xp, bz[j], ns[j])
                              
                misc.unpack(xp,bz,{'l':0,'q':[],'s':[ns[j]]},
                    offsetx = offsetj, offsety = offsetj)
               

                
                blas.tbmv(Gamma[j], bz, n = ns[j]**2, k = 0, 
                    ldA = 1, offsetx = offsetj) 
                
                # bz = Vt' * bz * Vt
                #    = uz

                cngrnc(Vt[j], bz, trans = 'T', offsetx = offsetj, n = ns[j])
                symmetrize(bz, ns[j], offset = offsetj)
                offsetj += ns[j]**2
                

            # x = -bz - r * uz * r' 
            blas.copy(z,x)
            blas.copy(bz,z)
            offsetj = 0
            for j in xrange(N):             
                cngrnc(+W['r'][j], bz, offsetx = offsetj, n = ns[j])
                offsetj += ns[j]**2
            blas.axpy(bz, x) 
            blas.scal(-1.0, x) 


        return solve
        
        
        
        
    sol = solvers.coneqp(Pf, q, Gf, h, dims, Af, b, kktsolver = F,  xdot=xdot)

    proxargs['X'][:] = +sol['s'][:]
    proxargs['dualy'][:] = +sol['y'][:]
    proxargs['dualS'][:] = +sol['z'][:]

    primal = blas.dot(proxargs['C'],proxargs['X'])
    gap = sol['gap']
    return primal, gap
 
def proxqp_clique_SNL(c, A, b, z, rho):

    """
    Solves the 1-norm regularized conic LP

        min.  < c, x > + || A(x) - b ||_1 + (rho/2) || x - z ||^2
        s.t.  x >= 0

    for a single dense clique .
    
    The method is used in this package to solve the 
    sensor node localization problem
    
    Input arguments.

        c is a 'd' matrix of size n_k**2 x 1

        A is a 'd' matrix.  with size n_k**2 times m_k.  Each of its columns 
            represents a symmetric matrix of order n_k in unpacked column-major 
            order. The term  A ( x ) in the primal constraint is given by

                A(x) = A' * vec(x).

            The adjoint A'( y ) in the dual constraint is given by

                A'(y) = mat( A * y ).

            Only the entries of A corresponding to lower-triangular
            positions are accessed.
             

        b is a 'd' matrix of size m_k x 1.

        z is a 'd' matrix of size n_k**2 x 1
        
        rho is a positive scalar.  


    Output arguments.

        sol : Solution dictionary for quadratic optimization problem.
        
        primal : objective for optimization problem without prox term (trace C*X)

    """
    
    ns2,ms = A.size
    nl,msl = len(b)*2,len(b)

    ns = int(sqrt(ns2))        
    dims = {'l': nl, 'q': [], 's': [ns] } 
    
    c = matrix([matrix(1.0,(nl,1)),c])
    z = matrix([matrix(0.0,(nl,1)),z])
    q = +c
    blas.axpy(z,q,alpha=-rho,offsetx=nl,offsety=nl)
    
    symmetrize(q,ns, offset = nl)
    q = q[:]
    h = matrix(0.0, (nl + ns2, 1))
    
    bz = +q
    xp = +q
    
    def P(u, v, alpha = 1.0, beta = 0.0):
        # v := alpha * rho * u + beta * v
        blas.scal(beta,v)
        blas.axpy(u,v,alpha=alpha*rho,offsetx = nl, offsety = nl)

    def xdot(x, y): 
        misc.trisc(x, dims)
        adot = blas.dot(x,y)
        misc.triusc(x, dims)
        return adot
     
    def Gf(u, v, alpha = 1.0, beta = 0.0, trans = 'N'):
        # v = -alpha*u + beta * v
        blas.scal(beta, v)
        blas.axpy(u,v,alpha=-alpha)
    
    def Af(u, v, alpha = 1.0, beta = 0.0,trans="N"):

        # v := alpha * A(u) + beta * v if trans is 'N'
        # v := alpha * A'(u) + beta * v if trans is 'T'
        blas.scal(beta,v)
        if trans == "N":
            blas.axpy(u,v,alpha=alpha,n=nl/2)
            blas.axpy(u,v,alpha=-alpha,offsetx = nl/2,n=nl/2)
            sgemv(A,u,v,n=ns,m=ms,alpha=alpha,beta=1.0,trans="T",offsetx = nl)
            
        elif trans == "T":
            blas.axpy(u,v,alpha=alpha,n=nl/2)
            blas.axpy(u,v,alpha=-alpha,offsety = nl/2,n=nl/2)
            sgemv(A,u,v,n=ns,m=ms,alpha=alpha,beta=1.0,trans="N",offsety = nl)
            
    
    U = matrix(0.0, (ns,ns)) 
    Vt = matrix(0.0, (ns,ns)) 
    sv = matrix(0.0, (ns,1)) 
    Gamma = matrix(0.0, (ns,ns))
    
    if type(A) is spmatrix:
        VecAIndex = +A[:].I
    As = matrix(A)
    Aspkd = matrix(0.0,((ns+1)*ns/2,ms))
    tmp = matrix(0.0,(ms,1))
    
    def F(W):
        # SVD R[j] = U[j] * diag(sig[j]) * Vt[j]
        lapack.gesvd(+W['r'][0], sv, jobu = 'A', jobvt = 'A',
            U = U, Vt = Vt)               
        
        
       
        W2 = mul(+W['d'],+W['d'])
               
        # Vt[j] := diag(sig[j])^-1 * Vt[j]
        for k in xrange(ns):
            blas.tbsv(sv, Vt, n = ns, k = 0, ldA = 1,offsetx = k*ns)

        # Gamma[j] is an ns[j] x ns[j] symmetric matrix
        #  (sig[j] * sig[j]') ./  sqrt(1 + rho * (sig[j] * sig[j]').^2)
        # S = sig[j] * sig[j]'
        S = matrix(0.0, (ns, ns))
        blas.syrk(sv, S)
        Gamma = div(S, sqrt(1.0 + rho * S**2))
        symmetrize(Gamma,ns)
  
        # As represents the scaled mapping
        #
        #     As(x) = A(u * (Gamma .* x) * u')
        #    As'(y) = Gamma .* (u' * A'(y) * u)
        #
        # stored in a similar format as A, except that we use packed 
        # storage for the columns of As[i][j].
         
        
        if type(A) is spmatrix:
            blas.scal(0.0,As)
            As[VecAIndex] = +A[VecAIndex]
        else:
            blas.copy(A,As)
            
                            
        # As[i][j][:,k] = diag( diag(Gamma[j]))*As[i][j][:,k]
        # As[i][j][l,:] = Gamma[j][l,l]*As[i][j][l,:]
        for k in xrange(ms):
            cngrnc(U, As, trans = 'T',offsetx=k*(ns2))
            blas.tbmv(Gamma, As, n = ns2,k = 0, ldA = 1,offsetx=k*(ns2))

        misc.pack(As,Aspkd,{'l': 0, 'q': [], 's': [ns]*ms})
        
        # H is an m times m block matrix with i, k block
        #
        #      Hik = sum_j As[i,j]' * As[k,j] 
        #
        # of size ms[i] x ms[k].  Hik = 0 if As[i,j] or As[k,j] 
        # are zero for all j
        H = matrix(0.0, (ms, ms)) 
        blas.syrk(Aspkd, H, trans = 'T',beta = 1.0, k = ns*(ns+1)/2)
        

        #H = H + spmatrix(W2[:nl/2] + W2[nl/2:] ,range(nl/2),range(nl/2))        
        blas.axpy(W2,H,n=ms,incy=ms+1,alpha=1.0)
        blas.axpy(W2,H,offsetx = ms,n=ms,incy=ms+1,alpha=1.0)

        lapack.potrf(H)
   
        
            
        def solve(x, y, z):
            
            """
            Returns solution of 

                rho * ux + A'(uy) - r^-T * uz * r^-1 = bx
                A(ux)                                = by
                -ux               - r * uz * r'      = bz.

            On entry, x = bx, y = by, z = bz.
            On exit, x = ux, y = uy, z = uz.
            """
            
            # bz is a copy of z in the format of x
            blas.copy(z,bz)
            blas.axpy(bz,x,alpha=rho,offsetx = nl, offsety = nl)
            # x := Gamma .* (u' * x * u)  
            #    = Gamma .* (u' * (bx + rho * bz) * u)

            cngrnc(U, x, trans = 'T',offsetx = nl)
            blas.tbmv(Gamma, x, n = ns2, k = 0, ldA = 1,offsetx=nl)
            blas.tbmv(+W['d'],x,n=nl,k=0,ldA=1)



            
            # y := y - As(x)
            #   := by - As( Gamma .* u' * (bx + rho * bz) * u)

            misc.pack(x, xp, dims)
            blas.gemv(Aspkd, xp, y, trans = 'T',alpha = -1.0, beta = 1.0, \
                m = ns*(ns+1)/2, n = ms,offsetx = nl)
            
            
            
            
            #y = y - mul(+W['d'][:nl/2],xp[:nl/2])+ mul(+W['d'][nl/2:nl],xp[nl/2:nl])
            blas.tbmv(+W['d'],xp,n=nl,k=0,ldA=1)
            blas.axpy(xp,y,alpha=-1,n=ms)
            blas.axpy(xp,y,alpha=1,n=ms,offsetx = nl/2)
        
            # y := -y - A(bz) 
            #    = -by - A(bz) + As(Gamma .*  (u' * (bx + rho * bz) * u) 
            

            Af(bz, y, alpha = -1.0, beta = -1.0)  

        
            
            
            # y := H^-1 * y
            #    = H^-1 ( -by - A(bz) + As(Gamma.* u'*(bx + rho*bz)*u) )
            #    = uy
            
            
          
            blas.trsv(H, y)
            blas.trsv(H, y,trans='T')


            # bz = Vt' * vz * Vt 
            #    = uz where
            # vz := Gamma .* ( As'(uy)  - x ) 
            #     = Gamma .* ( As'(uy)  - Gamma .* (u'*(bx + rho *bz)*u) )
            #     = Gamma.^2 .* ( u' * (A'(uy) - bx - rho * bz) * u ).
            
            misc.pack(x, xp, dims)
            blas.scal(-1.0, xp)
             
            
            
            blas.gemv(Aspkd, y, xp, alpha = 1.0,
                 beta = 1.0, m = ns*(ns+1)/2, n = ms,offsety=nl)

            #xp[:nl/2] = xp[:nl/2] + mul(+W['d'][:nl/2],y)
            #xp[nl/2:nl] = xp[nl/2:nl] - mul(+W['d'][nl/2:nl],y)
            
            
            blas.copy(y,tmp)
            blas.tbmv(+W['d'],tmp,n=nl/2,k=0,ldA=1)
            blas.axpy(tmp,xp,n=nl/2)

            blas.copy(y,tmp)
            blas.tbmv(+W['d'],tmp,n=nl/2,k=0,ldA=1,offsetA = nl/2)
            blas.axpy(tmp,xp,alpha=-1,n=nl/2,offsety=nl/2)
 
                
                 
            
            # bz[j] is xp unpacked and multiplied with Gamma
            blas.copy(xp,bz)#,n = nl)
            misc.unpack(xp, bz, dims)
            blas.tbmv(Gamma, bz, n = ns2, k = 0, ldA = 1,offsetx = nl) 

            # bz = Vt' * bz * Vt
            #    = uz
            cngrnc(Vt, bz, trans = 'T',offsetx = nl)
            
            symmetrize(bz,ns,offset = nl)
            

            # x = -bz - r * uz * r'
            # z contains r.h.s. bz;  copy to x
            #so far, z = bzc (untouched)
            blas.copy(z,x)
            blas.copy(bz, z)
                     
            cngrnc(W['r'][0], bz,offsetx=nl)
            blas.tbmv(W['d'], bz, n = nl, k = 0, ldA = 1) 

            blas.axpy(bz, x) 
            blas.scal(-1.0, x)

        return solve

    sol = solvers.coneqp(P, q, Gf, h, dims, Af, b, None, F, xdot=xdot)
    primal = blas.dot(sol['s'],c)

    sol['s'] = sol['s'][nl:]
    sol['z'] = sol['z'][nl:]

    return sol, primal


def proxqp_clique(c, A, b, z, rho):

    """
    Solves the conic QP

        min.  < c, x > + (rho/2) || x - z ||^2
        s.t.  A(x) = b
              x >= 0

    and its dual

        max.  -< b, y > - 1/(2*rho) * || c + A'(y) - rho * z - s ||^2 
        s.t.  s >= 0.

    for a single dense clique. 
    
    If the problem has block-arrow correlative sparsity, then the previous
    function 
    
    X = proxqp(c,A,b,z,rho,**kwargs)
    
    is equivalent to
    
    for k in xrange(ncliques):
        X[k] = proxqp_clique(c[k],A[k][k],b[k],z[k],rho,**kwargs)
    
    and each call can be implemented in parallel.
    
    Input arguments.

        c is a 'd' matrix of size n_k**2 x 1

        A is a 'd' matrix.  with size n_k**2 times m_k.  Each of its columns 
            represents a symmetric matrix of order n_k in unpacked column-major 
            order. The term  A ( x ) in the primal constraint is given by

                A(x) = A' * vec(x).

            The adjoint A'( y ) in the dual constraint is given by

                A'(y) = mat( A * y ).
             

        b is a 'd' matrix of size m_k x 1.
        
        z is a 'd' matrix of size n_k**2 x 1
        
        rho is a positive scalar.  

    Output arguments.
    
        sol : Solution dictionary for quadratic optimization problem.
        
        primal : objective for optimization problem without prox term (trace C*X)

    """

    ns2,ms = A.size
    ns = int(sqrt(ns2))
    dims = {'l':0,'q':[],'s':[ns]}
    
    q = +c
    blas.axpy(z,q,alpha=-rho)
    symmetrize(q,ns, offset = 0)
    q = q[:]
    h = matrix(0.0, (ns2, 1))
    

    bz = +q
    xp = +q

    def P(u, v, alpha = 1.0, beta = 0.0):
        # v := alpha * rho * u + beta * v
        #if not (beta==0.0):
        blas.scal(beta,v)
        blas.axpy(u,v,alpha=alpha*rho)

    def xdot(x, y): 
        misc.trisc(x, {'l': 0, 'q': [], 's': [ns]})
        adot = blas.dot(x,y)
        misc.triusc(x, {'l': 0, 'q': [], 's': [ns]})
        return adot
     
    def Gf(u, v, alpha = 1.0, beta = 0.0, trans = 'N'):

        # v = -alpha*u + beta * v
        # u and v are vectors representing N symmetric matrices in the
        # cvxopt format.
        blas.scal(beta, v)
        blas.axpy(u,v,alpha=-alpha)
    
    def Af(u, v, alpha = 1.0, beta = 0.0,trans="N"):

        # v := alpha * A(u) + beta * v if trans is 'N'
        # v := alpha * A'(u) + beta * v if trans is 'T'
        blas.scal(beta,v)
        if trans == "N":
            sgemv(A,u,v,n=ns,m=ms,alpha=alpha,beta=1.0,trans="T",offsetx = 0)
        elif trans == "T":
            sgemv(A,u,v,n=ns,m=ms,alpha=alpha,beta=1.0,trans="N",offsetx = 0)
            
    
    U = matrix(0.0, (ns,ns)) 
    Vt = matrix(0.0, (ns,ns)) 
    sv = matrix(0.0, (ns,1)) 
    Gamma = matrix(0.0, (ns,ns))
    
    
    if type(A) is spmatrix:
        VecAIndex = +A[:].I

    Aspkd = matrix(0.0,((ns+1)*ns/2,ms))
    As = matrix(A)

  
    
    def F(W):
        # SVD R[j] = U[j] * diag(sig[j]) * Vt[j]
        lapack.gesvd(+W['r'][0], sv, jobu = 'A', jobvt = 'A',
            U = U, Vt = Vt)               
        
       
               
        # Vt[j] := diag(sig[j])^-1 * Vt[j]
        for k in xrange(ns):
            blas.tbsv(sv, Vt, n = ns, k = 0, ldA = 1,offsetx = k*ns)

        # Gamma[j] is an ns[j] x ns[j] symmetric matrix
        #
        #  (sig[j] * sig[j]') ./  sqrt(1 + rho * (sig[j] * sig[j]').^2)

        # S = sig[j] * sig[j]'
        S = matrix(0.0, (ns, ns))
        blas.syrk(sv, S)
        Gamma = div(S, sqrt(1.0 + rho * S**2))
        symmetrize(Gamma,ns)
  
        # As represents the scaled mapping
        #
        #     As(x) = A(u * (Gamma .* x) * u')
        #    As'(y) = Gamma .* (u' * A'(y) * u)
        #
        # stored in a similar format as A, except that we use packed 
        # storage for the columns of As[i][j].
         

        
        if type(A) is spmatrix:            
            blas.scal(0.0,As)
            try:
                As[VecAIndex] = +A['s'][VecAIndex]
            except:
                As[VecAIndex] = +A[VecAIndex]
        else:
            blas.copy(A,As)     
                            
        # As[i][j][:,k] = diag( diag(Gamma[j]))*As[i][j][:,k]
        # As[i][j][l,:] = Gamma[j][l,l]*As[i][j][l,:]
        for k in xrange(ms):
            cngrnc(U, As, trans = 'T',offsetx=k*(ns2))
            blas.tbmv(Gamma, As, n = ns2,k = 0, ldA = 1,offsetx=k*(ns2))

        
        misc.pack(As,Aspkd,{'l': 0, 'q': [], 's': [ns]*ms})
        
        # H is an m times m block matrix with i, k block
        #
        #      Hik = sum_j As[i,j]' * As[k,j] 
        #
        # of size ms[i] x ms[k].  Hik = 0 if As[i,j] or As[k,j] 
        # are zero for all j
        H = matrix(0.0, (ms, ms)) 
        blas.syrk(Aspkd, H, trans = 'T',beta = 1.0, k = ns*(ns+1)/2)

       
       
        lapack.potrf(H)
   
        
            
        def solve(x, y, z):

            """
            Returns solution of 

                rho * ux + A'(uy) - r^-T * uz * r^-1 = bx
                A(ux)                                = by
                -ux               - r * uz * r'      = bz.

            On entry, x = bx, y = by, z = bz.
            On exit, x = ux, y = uy, z = uz.
            """
         
            # bz is a copy of z in the format of x
            blas.copy(z,bz)
            blas.axpy(bz,x,alpha=rho)

            # x := Gamma .* (u' * x * u)  
            #    = Gamma .* (u' * (bx + rho * bz) * u)

            cngrnc(U, x, trans = 'T',offsetx = 0)
            blas.tbmv(Gamma, x, n = ns2, k = 0, ldA = 1,offsetx=0)
            

            # y := y - As(x)
            #   := by - As( Gamma .* u' * (bx + rho * bz) * u)
            #blas.copy(x,xp)
            #pack_ip(xp,n = ns,m=1,nl=nl)
            misc.pack(x, xp, {'l': 0, 'q': [], 's': [ns]})
            
            blas.gemv(Aspkd, xp, y, trans = 'T',alpha = -1.0, beta = 1.0, \
                m = ns*(ns+1)/2, n = ms,offsetx = 0)
                 
                  
        
            # y := -y - A(bz) 
            #    = -by - A(bz) + As(Gamma .*  (u' * (bx + rho * bz) * u) 
            Af(bz, y, alpha = -1.0, beta = -1.0)  
            
        
            
            
            # y := H^-1 * y
            #    = H^-1 ( -by - A(bz) + As(Gamma.* u'*(bx + rho*bz)*u) )
            #    = uy
            
            blas.trsv(H, y)
            blas.trsv(H, y,trans='T')
          
            
            # bz = Vt' * vz * Vt 
            #    = uz where
            # vz := Gamma .* ( As'(uy)  - x ) 
            #     = Gamma .* ( As'(uy)  - Gamma .* (u'*(bx + rho *bz)*u) )
            #     = Gamma.^2 .* ( u' * (A'(uy) - bx - rho * bz) * u ).
            #blas.copy(x,xp)
            #pack_ip(xp,n=ns,m=1,nl=nl)
            
            misc.pack(x, xp, {'l': 0, 'q': [], 's': [ns]})
            blas.scal(-1.0, xp)

            
            
            blas.gemv(Aspkd, y, xp, alpha = 1.0,
                 beta = 1.0, m = ns*(ns+1)/2, n = ms,offsety=0)
               
            
            # bz[j] is xp unpacked and multiplied with Gamma
            misc.unpack(xp, bz, {'l': 0, 'q': [], 's': [ns]})
            blas.tbmv(Gamma, bz, n = ns2, k = 0, ldA = 1,offsetx = 0) 

            # bz = Vt' * bz * Vt
            #    = uz
            cngrnc(Vt, bz, trans = 'T',offsetx = 0)
            
            symmetrize(bz,ns,offset = 0)
         
            # x = -bz - r * uz * r'
            # z contains r.h.s. bz;  copy to x
            blas.copy(z,x)
            blas.copy(bz, z)
                     
            cngrnc(W['r'][0], bz,offsetx=0)
            blas.axpy(bz, x) 
            blas.scal(-1.0, x)
            
            
        return solve


    #solvers.options['show_progress'] = True
    sol = solvers.coneqp(P, q, Gf, h, dims, Af, b, None, F, xdot=xdot)
    primal = blas.dot(c,sol['s'])
    return sol, primal
    
    
def proxqp_blockdiagonal(proxargs,proxqp_per_clique):

    """
    Solves the conic QP

        min.  < c, x > + (rho/2) || x - z ||^2
        s.t.  A(x) = b
              x >= 0

    and its dual

        max.  -< b, y > - 1/(2*rho) * || c + A'(y) - rho * z - s ||^2 
        s.t.  s >= 0.
    
    where the constraints have block-diagonal correlative sparsity 
    (variables in block i only concern constraints in block i)

    The primal variable is x = matrix(x_0[:], ...,  x_{N-1}[:]) with x_k a symmetric
    matrix of order nk.
    The dual variables are y = matrix(y_0, ..., y_{M-1}), with y_k a column
    vector of length mk, and s = matrix(s_0[:], .., s_{N-1}[:]), with s_i a symmetric
    matrix of order nk.

    In the primal cost function, c = (c_0[:], ..., c_{N-1}[:]), with c_k a 
    symmetric matrix, and < c, x > = sum_j tr (c_j * x_j).
    In the dual cost function, b = matrix(b_0, ..., b_{M-1}), with b_k a 
    column vector, and < b, y > = sum_j b_i' * y_i.

    The mapping A(x) is defined as follows. 
    The value of u = A(x) is u = matrix(u_0, ..., u_{M-1}) with u_i a vector
    defined by

        u_i =  Ai(x_i),  i = 0, ..., M-1.

    The adjoint v = A'(y)  is v = (v_0, ..., v_{N-1}) with v_j a matrix
    defined as

        v_j = Aj'(y_j),  j = 0, ..., N-1.

    If we expand the primal and dual problems we therefore have

        min.  sum_j < c_j, x_j >  + (rho/2) * sum_j || x_j - z_j ||_F^2 
        s.t.  Ai(x_i) = b_i, i = 0, ..., M-1 
              x_j >= 0, j = 0, ..., N-1

    with variables x_j and its dual

        max.  sum_i b_i'*y_i - 1/(2*rho) *
              sum_j || c_j - Aj'(y_j) + rho *z_j - s_j ||_F^2
        s.t.  s_j >= 0, j = 0, ..., N-1

    with variables y_j, s_j.

    Input arguments.
    
        proxargs['C']   is a stacked vector containing vectorized 'd' matrices 
                        c_k of size n_k**2 x 1, representing symmetric matrices.

        proxargs['A']   is a list of 'd' matrices
            
            A[j] has size n_j**2 times m_i.  Each of its columns 
            represents a symmetric matrix of order n_j in unpacked column-major 
            order. The term  Ajj(x[j]) in the primal constraint is given by

                Ajj(x_j) = A[j]' * vec(x_j).

            The adjoint Ajj'(y[j]) in the dual constraint is given by

                Ajj'(y_j) = mat( [j] * y_j ).


        proxargs['b']   is a stacked vector containing constraint vectors of 
                        size m_i x 1.

        proxargs['Z']   is a stacked vector containing variable Z vectors of 
                        size n_k**2 x 1.
            
        proxargs['sigma'] is a positive scalar (step size).
        
        proxargs['X']   contains the primal variable X in stacked vector form. 
        
        proxargs['dualy']   contains the dual variable y 
        
        proxargs['dualS']   contains the primal variable S
	
	    On output, proxargs['X'], proxargs['dualy'], and proxargs['dualZ'] 
	    will be filled with prox optimal primal and dual variables.
	    
        
    kwargs: 
        
        ns: List containing the order of each semidefinite subproblem
        
        ms : list of constrant block sizes
            
        multiprocess : Number of processors available on machine. In 

            multiprocessing mode, each clique is round-robin assigned to each 
            processor, and each prox is evaluated independently.
            
            If multiprocess = 1, a serial implementation is used.
	    
    Output arguments.

        primal : trace(C*X) where X is the optimal variable
        
        gap : trace(X*S) as computed by  CVXOPT
	

  
        
    OUTPUT:
    
        F : function pointer for prox evaluations
  
    """
    
    C, A, b = proxargs['C'],proxargs['A'],proxargs['b']
    ns, ms = proxargs['ns'],proxargs['ms']  
    multiprocess = proxargs['multiprocess']

    
    
    def F(proxargs):
        """
        Solves the conic QP

            min.  < c, x > + (rho/2) || x - z ||^2
            s.t.  A(x) = b
                  x >= 0

        and its dual

            max.  -< b, y > - 1/(2*rho) * || c + A'(y) - rho * z - s ||^2 
            s.t.  s >= 0.
            
        INPUT:
            
            proxargs : dictionary containing problem parameters and variables
        
        OUTPUT: 
            
            primal : Primal optimal value (trace (C*X))
            
            gap : Duality gap evaluated by CVXOPT (trace(X*S))
            
        """
        
        def proxqp_multiprocess(C,A,b,Z,rho,q):
            obj,pv = proxqp_per_clique(C,A,b,Z, rho)
            q.put([+obj['s'],+obj['y'],+obj['z'], pv,obj['gap']])
        rho = proxargs['sigma'] 
        Z, X = proxargs['Z'],proxargs['X']
        y, S = proxargs['dualy'], proxargs['dualS']
        primal,gap = 0.,0.
        if multiprocess == 1:
            offsetns,offsetms = 0,0

            for k in xrange(len(ns)):            
                obj,pv = proxqp_per_clique(C[offsetns : offsetns + ns[k]**2],A[k],
                    b[offsetms: offsetms + ms[k]], Z[offsetns : offsetns + ns[k]**2], 
                    rho)
                gap += obj['gap']
                primal += pv
                X[offsetns : offsetns + ns[k]**2] = +obj['s']
                y[offsetms:(offsetms+ms[k])] = +obj['y']
                S[offsetns : offsetns + ns[k]**2] = +obj['z']
                offsetms += ms[k]
                offsetns += ns[k]**2
        else:
            offsetns,offsetms = 0,0
            proclist = [[] for k in xrange(multiprocess)]
            nslist, mslist = [],[]
            for k in xrange(len(ns)):           
                q = Queue()
                proc = Process(target = proxqp_multiprocess, \
                        args=(C[offsetns : offsetns + ns[k]**2], A[k],  \
                        b[offsetms: offsetms + ms[k]], \
                        Z[offsetns : offsetns + ns[k]**2], rho,q))
                proc.start()
                proclist[k % multiprocess].append((k,proc,q))
                nslist.append(offsetns)
                mslist.append(offsetms)
                offsetms += ms[k]
                offsetns += ns[k]**2
            while(1):
                alldone = True
                for i in xrange(multiprocess):
                    if len(proclist[i]) > 0:
                        alldone = False
                        k = proclist[i][0][0]
                        xys = proclist[i][0][2].get()
                        xx,yy,ss,pp,gg = xys[0],xys[1],xys[2],xys[3],xys[4]
                        X[nslist[k]:(nslist[k]+ns[k]**2)] = xx
                        y[mslist[k]:(mslist[k]+ms[k])] = yy
                        S[nslist[k]:(nslist[k]+ns[k]**2)] = ss
                        gap += gg
                        primal += pp
                        proclist[i][0][1].join()
                        proclist[i].pop(0)
                if alldone: break
        return primal , gap
    return F

