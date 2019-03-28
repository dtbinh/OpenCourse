from cvxopt import spmatrix, sparse, matrix, blas, lapack, misc, base

def tril(A):
    """
    Returns the lower triangular part of A.
    """
    if isinstance(A,spmatrix):
        idx = [i for i,ij in enumerate(zip(A.I,A.J)) if ij[0]>=ij[1]]
        return spmatrix(A.V[idx],A.I[idx],A.J[idx],A.size)
    elif isinstance(A,matrix):
        B = matrix(0.0, A.size)
        lapack.lacpy(A, B, uplo = 'L')
        return B
    else:
        raise TypeError

def triu(A):
    """
    Returns the upper triangular part of A.
    """
    if isinstance(A,spmatrix):
        idx = [i for i,ij in enumerate(zip(A.I,A.J)) if ij[0]<=ij[1]]
        return spmatrix(A.V[idx],A.I[idx],A.J[idx],A.size)
    elif isinstance(A,matrix):
        B = matrix(0.0, A.size)
        lapack.lacpy(A, B, uplo = 'U')
        return B
    else:
        raise TypeError

def symmetrize(x,n = None,offset = 0):
    """
    Copies lower triangular part to upper triangular part
    """ 
    if n is None:  n = x.size[0]
    if type(x) is spmatrix:
        idx = [i for i,ij in enumerate(zip(x.I,x.J)) if ij[0] > ij[1]]
        return tril(x) + spmatrix(x.V[idx], x.J[idx], x.I[idx], x.size)
    else:
        for j in range(n-1):  
            blas.copy(x,x,offsetx = j*n+j+1+offset,offsety = j*n+j+n+offset,n=n-1-j,incy = n)
        


def perm(A, p):
    """
    Symmetric permutation of a symmetric sparse matrix. 

    :param A:    spmatrix
    :param p:    permutation vector

    """
    assert A.size[0] == A.size[1], "A must be a square matrix"
    assert A.size[0] == len(p), "length of p must be equal to the order of A"
    return A[p,p]

def eye(n):
    """
    Identity matrix of order n.
    """
    I = matrix(0.0,(n,n))
    I[::n+1] = 1.0
    return I
    
def reshape(X,N,M):
    assert (X.size[0]*X.size[1]) == (N*M), "reshape requires matrix dimension to stay the same!"
    if type(X) is matrix: return matrix(X,(N,M))
    X = sparse(X[:])
    I,V = X.I,X.V
    J = [(i / N) for i in list(I)]
    I = [(i % N) for i in list(I)]
    return spmatrix(V,I,J,(N,M))
def spdiag_rect(Xlist):
    """
    Creates a sparse block-diagonal matrix where each block may be non-square.

    ARGUMENTS
        Xlist : list of blocks
    
    RETURNS
        X : sparse block-diagonal matrix
    """
    I,J,V = [],[],[]
    offsetx, offsety = 0,0
    
    for x in Xlist:
        x = sparse(x)
        I.append(x.I+offsety)
        J.append(x.J+offsetx)
        V.append(x.V)
        offsety += x.size[0]
        offsetx += x.size[1]
    I = matrix(I)
    J = matrix(J)
    V = matrix(V)
    return spmatrix(V,I,J,size=(offsety,offsetx))
    

def cngrnc(r, x, alpha = 1.0, trans = 'N', offsetx = 0, n = None):
    """
    In-place congruence transformation

        x := alpha * r * x * r' if trans = 'N'.

        x := alpha * r' * x * r if trans = 'T'.

    r is a square n x n-matrix. 
    x is a square n x n-matrix or n^2-vector. 
    """

    if n is None: n = r.size[0]

    # Scale diagonal of x by 1/2.  
    blas.scal(0.5, x, n = n, inc = n+1, offset = offsetx)
    
    # a := r*tril(x) if trans is 'N',  a := tril(x)*r if trans is 'T'.
    a = +r
    if trans == 'N':  
        blas.trmm(x, a, side = 'R', m = n, n = n, ldA = n, ldB = n, 
            offsetA = offsetx)
    else:  
        blas.trmm(x, a, side = 'L', m = n, n = n, ldA = n, ldB = n,
            offsetA = offsetx)

    # x := alpha * (a * r' + r * a')  
    #    = alpha * r * (tril(x) + tril(x)') * r' if trans is 'N'.
    # x := alpha * (a' * r  + r' * a)  
    #    = alpha * r' * (tril(x)' + tril(x)) * r if trans = 'T'.
    blas.syr2k(r, a, x, trans = trans, alpha = alpha, n = n, k = n, 
        ldB = n, ldC = n, offsetC = offsetx)


def sgemv(A, x, y, trans = 'N', alpha = 1.0, beta = 0.0, n = None, 
    m = None, offsetx = 0, offsety = 0):
    """
    A is an n^2 x m matrix with columns interpreted as vectorized
    symmetric n x n matrices in 'L' storage.

    If trans is 'N':

        y := alpha * A * x + beta * y

    y is a vector of length n**2 or an nxn-matrix representing a symmetric
    matrix in 'L' storage.  x is a vector of length m.

    If trans is 'T':
 
        y := alpha * A' * x + beta * y

    y is a vector of length m.   x is a vector of length n**2 or an 
    nxn-matrix representing an n x n symmetric matrix in 'L' storage

    The default values of the dimensions are n = int(sqrt(A.size[0])) and
    m = A.size[1].
    """

    if m is None: m = A.size[1]
    if n is None: n = int(sqrt(A.size[0]))
    if trans == 'T' and alpha:  misc.trisc(x, {'l': offsetx, 'q': [], 's': [n]})
    if type(A) is spmatrix:
        base.gemv(A, x, y, trans = trans, alpha = alpha, beta = beta, m = n**2,
            n = m,offsetx = offsetx, offsety = offsety) 
    else:
        blas.gemv(A, x, y, trans = trans, alpha = alpha, beta = beta, m = n**2,
            n = m,offsetx = offsetx, offsety = offsety) 
    if trans == 'T' and alpha: misc.triusc(x, {'l': offsetx, 'q': [], 's': [n]})


