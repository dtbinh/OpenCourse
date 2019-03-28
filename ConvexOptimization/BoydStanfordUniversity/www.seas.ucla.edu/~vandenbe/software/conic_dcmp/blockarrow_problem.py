from cvxopt import matrix,solvers, spmatrix, setseed, normal, lapack
from misc import eye

def blockarrow_problem(p,r,s,arrowwidth,seed=None):

    """
    Generates problem variables A, b, C, L for converted and unconverted
    problems. In the unconverted form, the aggregate sprasity is
    block-diagonal-arrow, and in the converted form, the correlative 
    sparsity is block-diagonal.
    
    ARGUMENTS:
        p : order of each diagonal block
        
        r : number of diagonal blocks
        
        s : number of constraints per clique
        
        arrowwidth : width of arrow in aggregate sparsity pattern
        
        seed : random generator seed
    
    RETURNS: dictionary with following fields
    
        Au, bu, Cu, X0u : unconverted problem parameters and feasible starting
            point X0u. Each is given as a matrix.
                    
        Sp : Aggregate sparsity pattern, represented by binary matrix
    
    """
    def get_blockarrow_index(k,p,r,arrowwidth):
        """
        Macro for getting indices of kth blockarrow clique.
        """
        index = range(k*p,(k+1)*p)
        index.extend(range(r*p,r*p+arrowwidth))
        return index
    
    def get_psd_matrix(p):
        tmp = matrix(normal((p)**2),(p,p))/2.0
        tmp = tmp + tmp.T
        while(1):
            try:
                lapack.potrf(+tmp)
                break
            except:
                tmp = tmp + .1*eye(p)
        return tmp
        
    setseed(seed)
    
    NVar = p*r+arrowwidth
    NCon = s*r
    
    # Generate sparsity pattern
    
    I = []
    J = []
    val = []

    for k in xrange(r):
        I.extend(range(p*k,p*k+p)*(p))
        ji = [[i]*(p) for i in range(p*k,p*k+p)]
        J.extend([i for sublist in ji for i in sublist])
    for k in xrange(r):
        ji = [[i]*(arrowwidth) for i in range(p*k,p*k+p)]
        I.extend(range(r*p,r*p+arrowwidth)*(p))
        J.extend([i for sublist in ji for i in sublist])
        J.extend(range(r*p,r*p+arrowwidth)*(p))
        I.extend([i for sublist in ji for i in sublist])
    k = r
    ji = [[i]*(arrowwidth) for i in range(p*k,p*k+arrowwidth)]
    I.extend(range(r*p,r*p+arrowwidth)*(arrowwidth))
    J.extend([i for sublist in ji for i in sublist])
    
    # Generate X and S
    
    X = spmatrix(0.0,I,J,size=(NVar,NVar))
    S = spmatrix(0.0,I,J,size=(NVar,NVar))
    
    for k in xrange(r):
        index = get_blockarrow_index(k,p,r,arrowwidth)
        tmp = get_psd_matrix(p+arrowwidth)
        X[index,index] = tmp
        tmp = get_psd_matrix(p+arrowwidth)
        S[index,index] = S[index,index] + tmp
    Sp = spmatrix(1,I,J,size=(NVar,NVar))    
    # Generate A and b
    b = matrix(0.0,(NCon,1));
    AI,AJ,AV = [],[],[]
    C = +S
    for k in xrange(r):
        index = get_blockarrow_index(k,p,r,arrowwidth)
        block = []
        for i in xrange(s):
            square = matrix(normal((p+arrowwidth)**2),(p+arrowwidth,p+arrowwidth))/2.0
            square = square + square.T
            C[index,index] = C[index,index] + normal(1)[0]*square
            block.append(square[:].T)
        block = matrix(block).T   
        b[s*k:(k+1)*s] =  block.T*X[index,index][:]
        for i in xrange(s):
            AI.extend([(l+NVar*j)  for j in index for l in index])
            AJ.extend([(s*k+i)] * ((p+arrowwidth)**2))
        AV.extend(list(block[:]))


    A = spmatrix(AV,AI,AJ,size=(NVar**2,NCon))
    return A,b,C,Sp

