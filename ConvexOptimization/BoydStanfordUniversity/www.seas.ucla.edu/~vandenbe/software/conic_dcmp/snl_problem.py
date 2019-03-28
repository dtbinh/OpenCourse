from cvxopt import spmatrix, blas, matrix, sparse, amd, normal, lapack, sqrt,misc,  div, mul, spdiag
import random



def get_EDM(X,neighbors):
    """
    Finds sparse EDM 
    
    PARAMS:
        X : N x k sized data matrix 
            N = number of sensors, 
            k = data dimension.
        neighbors : number of nearest neighbors
        
    RETURNS
        EDM : partial Euclidean distance matrix where
        
            EDM[i,j] = ||Xi-Xj||_F^2 if i is a nearest neighbor of j OR 
                    j is a nearest neighbor of i
            EDM[i,j] = 0 otherwise
            
        NN : adjacency matrix for nearest neighbor connections
    """
    n,k = X.size
    dG = mul(X,X)*matrix(1.,(k,1))
    I,J,V = [],[],[]
    if n <= 100:
        G = X*X.T
        o1 = matrix(list(dG)*n,(n,n))
        o2 = [[dG[i]]*n for i in xrange(n)]
        o2 = matrix([i for j in o2 for i in j],(n,n))
        EDM = o1+o2-2*G
        for i in xrange(n):
            tmp = list(EDM[:,i])
            tmp.sort()
            
            thresh = tmp[neighbors + 1]
            for j in xrange(n):
                if EDM[j,i] <= thresh: 
                    I.append(j)
                    J.append(i)
                    V.append(EDM[j,i])
    else:
        block = 100
        o1 = matrix(list(dG)*block,(n,block))
        offset = 0
        while(offset < n):
            if ((offset + block) > n): block = n - offset
            o2 = [[dG[i]]*n for i in xrange(offset,offset+block)]
            o2 = matrix([i for j in o2 for i in j],(n,block))
            G = X*X[offset : (offset + block),:].T

            EDMsect = o1 + o2 - 2*G
            
            for i in xrange(block):
                tmp = list(EDMsect[:,i])
                tmp.sort()
                
                thresh = tmp[neighbors + 1]
                for j in xrange(n):
                    if EDMsect[j,i] <= thresh: 
                        I.append(j)
                        J.append(i + offset)
                        V.append(EDMsect[j,i])
        

            offset += block
    EDM = spmatrix(V,I,J,(n,n))    
    for k in xrange(len(I)):
        i = I[k]
        j = J[k]
        if EDM[i,j] > 0: EDM[j,i] = EDM[i,j]
        if EDM[j,i] > 0: EDM[i,j] = EDM[j,i]
    NN = spmatrix(1,EDM.I,EDM.J,EDM.size)
    return EDM, NN
    

def snl_problem(seed=None, **kwargs):

    """
    Uniformly scatters N sensors in unit hypercube and generates partial EDM matrix
    
    ARGUMENTS:
        
        seed : random generator seed
        
        kwargs:
            n : number of sensors
            
            k : dimension of space
            
            neighbors : number of nearest neighbors
    
            noise : measurement distance noise
            
            
    RETURNS: 
    
        X : n x k data matrix where X[i,:] is the position of sensor i
        
        EDM : partial Euclidean distance matrix where nearest neighbor distances
                are filled
                
        NN : adjacency matrix for nearest neighbor distances
    
    """
    
    n = kwargs.get('nsensors',100) 
    k = kwargs.get('dim',3)
    neighbors = kwargs.get('neighbors',5)
    noise = kwargs.get('noise',0.0)
    
    # generating problem data
    random.seed(seed)
    X = matrix([random.uniform(0,1) for i in xrange(n*k)],(n,k))
    EDM,NN = get_EDM(X,neighbors)
    
    I,J,V = EDM.I, EDM.J, EDM.V
    EDMnoise = spmatrix(normal(len(I),1),I,J,(n,n))
    EDMnoise = (EDMnoise + EDMnoise.T)/2
    EDM[:] = mul(EDM,(1.0+ noise*EDMnoise))[:]
    for i in xrange(n): EDM[i,i] = 0.0
    
    return X, EDM, NN
    
def get_unconverted_problem(X, EDM, Spun):
    """
    Generates problem variables A, b, C, L for converted and unconverted
    problems. The aggregate sparsity pattern of the unconverted problem

    is equivalent to the nearest neighbor adjacency matrix.
    In the converted form, the correlative 
    sparsity is block-diagonal.
    
    ARGUMENTS:

        X : n x k data matrix where X[i,:] is the position of sensor i
        
        EDM : partial Euclidean distance matrix where nearest neighbor distances
                are filled
                
        NN : adjacency matrix for nearest neighbor distances
    
            
    RETURNS: 
    
        Aun, bun, Cun : unconverted problem parameters.
        
        dims : dictionary of cone dimensions, where
            
            dims['l'] : the number of nonnegative linear cones (= 2*number of edges)
            
            dims['q'] : list of quadratic cone orders (none)
            
            dims['s'] : list of semidefinite cone orders (only one element = n)   
    
    """
    
    n = EDM.size[0]

    
    # get unconverted
    Aun = ([],[],[],[],[],[])
    bun = []

    nzsquare = []

    ncol = 0
    kl = 0
    for j in xrange(n):
        kl += j
        for i in xrange(j+1,n):
            if Spun[i,j] > 0:
                k = i + n*j
                Aun[0].extend([1.0,1.0,-1.0,-1.0])
                Aun[1].extend([i + n*i,j + n*j,i + n*j,j + n*i])
                Aun[2].extend([ncol,ncol,ncol,ncol])
                Aun[3].extend([1.,-1.])
                Aun[4].extend([ncol])
                Aun[5].extend([ncol,ncol])
                ncol += 1
                bun.append(EDM[i,j])
    I2 = []
    for k in Aun[4]:
        I2.extend([Aun[4][k],Aun[4][k]+ncol])
    Aunl = spmatrix(Aun[3],I2,Aun[5],size=(ncol*2,ncol))
    Auns = spmatrix(Aun[0],Aun[1],Aun[2],size=(n**2,ncol))
    Aun = sparse([Aunl,Auns])
    
    bun = matrix(bun)

    Cun = matrix([matrix(1.,(2*ncol,1)),matrix(0.0,(n**2,1))])
    dims = {'l':ncol*2,'q':[],'s':[n]}
    return Aun,bun,Cun,dims
