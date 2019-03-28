from cvxopt import matrix, spmatrix, normal, blas, printing
from misc import tril, perm, symmetrize
from types import BuiltinFunctionType, FunctionType

def lmerge(left, right, offsetl, offsetr, nl, nr):
    
    tmp = matrix(0,(nl+nr,1))
    il = 0
    ir = 0
    k = 0
    k1 = k;
    while (il < nl and ir < nr):
        if (left[offsetl+il] < right[offsetr+ir]):
            tmp[k] = left[offsetl+il];
            il += 1
        elif (left[offsetl+il] > right[offsetr+ir]):
          tmp[k] = right[offsetr+ir]
          ir += 1
        else:
            tmp[k] = left[offsetl+il]
            il += 1
            ir += 1
        k += 1
    if (il < nl) :
        for i in xrange(nl-il): tmp[k+i] = left[offsetl+il+i]
        k += nl-il

    if (ir < nr):
        for i in xrange(nr-ir): tmp[k+i] = right[offsetr+ir+i]
        k += nr-ir
    for i in xrange(k): 
        left[offsetl+i] = tmp[i]
    
    return k
    

def __tdfs(j, k, head, next, post, stack):
    """
    Depth-first search and postorder of a tree rooted at node j.
    """
    top = 0
    stack[0] = j
    while (top >= 0):
        p = stack[top]
        i = head[p]
        if i == -1:
            top -= 1
            post[k] = p
            k += 1
        else:
            head[p] = next[i]
            top += 1
            stack[top] = i

    return k

def post_order(parent):
    """
    Post order a forest.
    """
    n = len(parent)
    k = 0

    p = matrix(0,(n,1))
    head = matrix(-1,(n,1))
    next = matrix(0,(n,1))
    stack = matrix(0,(n,1))

    for j in range(n-1,-1,-1):
        if (parent[j] == -1): continue
        next[j] = head[parent[j]]
        head[parent[j]] = j

    for j in range(n):
        if (not parent[j] == -1): continue
        k = __tdfs(j, k, head, next, p, stack)

    return p

def etree(A):
    """
    Compute elimination tree from upper triangle of A.
    """
    assert isinstance(A,spmatrix), "A must be a sparse matrix"
    assert A.size[0] == A.size[1], "A must be a square matrix"

    n = A.size[0]
    cp,ri,_ = A.CCS
    parent = matrix(0,(n,1))
    w = matrix(0,(n,1))

    for k in range(n):
        parent[k] = -1
        w[k] = -1
        for p in range(cp[k],cp[k+1]):
            i = ri[p]
            while ((not i == -1) and (i < k)):
                inext = w[i]
                w[i] = k
                if inext == -1: parent[i] = k
                i = inext;

    return parent

def __leaf(i, j, first, maxfirst, prevleaf, ancestor):
    """
    Determine if j is leaf of i'th row subtree.
    """
    jleaf = 0
    if i<=j or first[j] <= maxfirst[i]: return -1, jleaf
    maxfirst[i] = first[j]
    jprev = prevleaf[i]
    prevleaf[i] = j
    if jprev == -1: jleaf = 1
    else: jleaf = 2
    if jleaf == 1: return i, jleaf
    q = jprev
    while q != ancestor[q]: q = ancestor[q]
    s = jprev
    while s != q:
        sparent = ancestor[s]
        ancestor[s] = q
        s = sparent
    return q, jleaf

def counts(A, parent, post):
    """
    Compute column counts.
    """
    n = A.size[0]
    colcount = matrix(0,(n,1))
    ancestor = matrix(range(n),(n,1))
    maxfirst = matrix(-1,(n,1))
    prevleaf = matrix(-1,(n,1))
    first = matrix(-1,(n,1))
    
    for k in range(n):
        j = post[k]
        if first[j] == -1:
            colcount[j] = 1
        else:
            colcount[j] = 0
        while  j != -1 and first[j] == -1:
            first[j] = k;
            j = parent[j]

    cp,ri,_ = A.CCS
    for k in range(n):
        j = post[k]
        if parent[j] != -1:
            colcount[parent[j]] -= 1
        for p in range(cp[j],cp[j+1]):
            i = ri[p]
            if i <= j: continue
            q, jleaf = __leaf(i, j, first, maxfirst, prevleaf, ancestor)
            if jleaf >= 1: colcount[j] += 1
            if jleaf == 2: colcount[q] -= 1
        if parent[j] != -1: ancestor[j] = parent[j]
    for j in range(n):
        if parent[j] != -1: colcount[parent[j]] += colcount[j]

    return colcount

def pothen_sun(par, post, colcount):
    """
    Find supernodes and supernodal etree.

    ARGUMENTS
    par       parent array

    post      array with post ordering

    colcount  array with column counts

    RETURNS
    snpar     supernodal parent structure 

    flag      integer vector of length n; if flag[i] < 0, then -flag[i]
              is the degree of the supernode with repr. vertex i; if
              flag[i] >= 0, then flag[i] is the repr. vertex to which
              node i belongs.
    """
    
    n = len(par)
    flag = matrix(-1, (n, 1))
    snpar = matrix(-1, (n, 1))
    snodes = n
    ch = {}
    
    for j in post:

        if par[j] in ch: ch[par[j]].append(j)
        else: ch[par[j]] = [j]

        mdeg = colcount[j] - 1

        if not par[j] == -1:
            if mdeg == colcount[par[j]] and flag[par[j]] == -1:
                # par[j] not assigned to supernode
                snodes -= 1
                if flag[j] < 0:   # j is a repr. vertex
                    flag[par[j]] = j
                    flag[j] -= 1
                else:             # j is not a repr. vertex
                    flag[par[j]] = flag[j]
                    flag[flag[j]] -= 1
        else:
            if flag[j] < 0: snpar[j] = -1
            else: snpar[flag[j]] = -1

        if flag[j] < 0: k = j
        else: k = flag[j]

        if j in ch:
            for i in ch[j]:
                if flag[i] < 0: l = i
                else: l = flag[i]
                if not l == k: snpar[l] = k


    repr = matrix([i for i in range(n) if flag[i] < 0])
    deg = matrix([-flag[i] for i in range(n) if flag[i] < 0])

    # renumber etree with number of supernodes
    sn = matrix(-1, (n+1, 1))
    for k, r in enumerate(repr): sn[r] = k
    snpar = sn[snpar[repr]]

    return snpar, flag

def supernodes(par, post, colcount):
    """
    Find supernodes.

    ARGUMENTS
    par       parent array

    post      array with post ordering

    colcount  array with column counts

    RETURNS
    snode     array with supernodes; snode[snptr[k]:snptr[k+1]] contains
              the indices of supernode k

    snptr     pointer array; snptr[k] is the index of the representative
              vertex of supernode k in the snode array

    snpar     supernodal parent structure 
    """
    snpar, flag = pothen_sun(par, post, colcount)
    n = len(par)
    N = len(snpar)

    snode = matrix(0, (n,1))
    snptr = matrix(0, (N+1,1))

    slist = [[] for i in range(n)]
    for i in range(n):
        f = flag[i]
        if f < 0:
            slist[i].append(i)
        else:
            slist[f].append(i)

    k = 0; j = 0
    for i,sl in enumerate(slist):
        nsl = len(sl)
        if nsl > 0:
            snode[k:k+nsl] = matrix(sl)
            snptr[j+1] = snptr[j] + nsl
            k += nsl
            j += 1
        
    return snode, snptr, snpar

def amalgamate(colcount, snode, snptr, snpar, snpost, merge_function):
    """
    Supernodal amalgamation.

       colcount, snode, snptr, snpar, snpost = ...
         amalgamate(colcount, snode, snptr, snpar, snpost, merge_function)
    
    PURPOSE
    Iterates over the clique tree in topological order and greedily
    merges a supernode with its parent if

       merge_function(|J_{par(k)}|, |J_k|, |N_{par(k)}|, |N_k|)

    returns True.

    ARGUMENTS
    colcount  vector with column counts

    snode     vector with supernodes
 
    snptr     vector with offsets

    snpar     vector with supernodal parent indices

    snpost    vector with supernodal post ordering

    merge_function
              function 

    RETURNS
    colcount  vector with amalgamated column counts

    snode     vector with amalgamated supernodes
 
    snptr     vector with amalgamated offsets

    snpar     vector with amalgamated supernodal parent indices

    snpost    vector with amalgamated supernodal post ordering
    """
    N = len(snpost)
    ch = {}
    for j in snpost:
        if snpar[j] in ch: ch[snpar[j]].append(j)
        else: ch[snpar[j]] = [j]

    snlist = []
    for k in range(N):
        snlist.append(snode[snptr[k]:snptr[k+1]])

    snpar_ = +snpar
    colcount_ = +colcount
    Ns = N
    for k in snpost:
        if snpar_[k] != -1:
            colk = colcount_[snlist[k][0]]
            colp = colcount_[snlist[snpar_[k]][0]]
            nk = len(snlist[k])
            np = len(snlist[snpar_[k]])
            if merge_function and merge_function(colp,colk,np,nk):
                # merge supernode k and snpar[k]
                snlist[snpar_[k]] = matrix(sorted(list(snlist[k]) + list(snlist[snpar_[k]])))
                snlist[k] = None
                colcount_[snlist[snpar_[k]][0]] = colp + nk
                Ns -= 1
                if ch.has_key(k):
                    for c in ch[k]:
                        snpar_[c] = snpar_[k]
                    ch[snpar_[k]] += ch[k]
                snpar_[k] = -1

    L = [i for i,s in enumerate(snlist) if s is not None]
    snptr_ = +snptr
    snode_ = +snode
    for i,l in enumerate(L):
        snptr_[i+1] = snptr_[i] + len(snlist[l])
        snode_[snptr_[i]:snptr_[i+1]] = snlist[l]

    snpar_ = snpar_[L]
    for i in range(len(snpar_)):
        if snpar_[i] != -1:
            snpar_[i] = L.index(snpar_[i])
    snpost_ = post_order(snpar_)
    return colcount_, snode_, snptr_, snpar_, snpost_

def embed(A, colcount, snode, snptr, snpar, snpost):
    """
    Compute filled pattern.

       colptr, rowidx = embed(A, colcount, snode, snptr, snpar, snpost)

    PURPOSE
    Computes rowindices and column pointer for representative vertices in supernodes.

    ARGUMENTS
    A         sparse matrix

    colcount  vector with column counts

    snode     vector with supernodes
 
    snptr     vector with offsets

    snpar     vector with supernodal parent indices

    snpost    vector with supernodal post ordering

    RETURNS
    colptr    vector with offsets 

    rowidx    vector with rowindices 
    """

    Alo = tril(A)
    cp,ri,_ = Alo.CCS
    N = len(snpar)

    # colptr for compressed cholesky factor
    colptr = matrix(0,(N+1,1))
    for k in range(N):
        colptr[k+1] = colptr[k] + colcount[snode[snptr[k]]]
    rowidx = matrix(-1,(colptr[-1],1))
    cnnz = matrix(0,(N,1))

    # compute compressed sparse representation
    for k in range(N):
        p = snptr[k]
        Nk = snptr[k+1]-p
        nk = cp[snode[p]+1] - cp[snode[p]]
        rowidx[colptr[k]:colptr[k]+nk] = ri[cp[snode[p]]:cp[snode[p]+1]]
        cnnz[k] = nk
        for i in range(1,Nk):
            nk = cp[snode[p+i]+1]-cp[snode[p+i]]
            cnnz[k] = lmerge(rowidx, ri, colptr[k], cp[snode[p+i]], cnnz[k], nk)

    for k in snpost:
        p = snptr[k]
        Nk = snptr[k+1]-p
        if snpar[k] != -1:
            cnnz[snpar[k]] = lmerge(rowidx,rowidx,colptr[snpar[k]], colptr[k]+Nk,cnnz[snpar[k]], cnnz[k]-Nk)

    return colptr, rowidx

def relative_idx(colptr, rowidx, snptr, snpar):
    """
    Compute relative indices of update matrices in frontal matrix of parent.
    """
    
    relptr = matrix(0, (len(snptr),1))
    relidx = matrix(-1, (colptr[-1],1))

    def lfind(a,b):
        i = 0
        ret = +a
        for k in range(len(a)):
            while a[k] != b[i]: i += 1
            ret[k] = i
            i += 1
        return ret
    
    for k in range(len(snpar)):
        p = snpar[k]
        relptr[k+1] = relptr[k]
        if p != -1:
            nk = snptr[k+1] - snptr[k]
            relptr[k+1] += colptr[k+1] - colptr[k] - nk
            relidx[relptr[k]:relptr[k+1]] = lfind(rowidx[colptr[k]+nk:colptr[k+1]], rowidx[colptr[p]:colptr[p+1]])

    return relptr, relidx[:relptr[k+1]]

def peo(A, p):
    """
    Checks whether an ordering is a perfect elmimination order.

    Returns `True` if the permutation :math:`p` is a perfect elimination order
    for a Cholesky factorization :math:`PAP^T = LL^T`. Only the lower
    triangular part of :math:`A` is accessed.
    
    :param A:   :py:class:`spmatrix`
    
    :param p:   :py:class:`matrix` or :class:`list` of length `A.size[0]`
    """

    n = A.size[0]
    assert type(A) == spmatrix, "A must be a sparse matrix"
    assert A.size[1] == n, "A must be a square matrix"
    assert len(p) == n, "length of p must be equal to the order of A" 

    As = symmetrize(A)
    cp,ri,_ = As.CCS

    # compute inverse permutation array
    ip = matrix(0,(n,1))
    ip[p] = matrix(range(n),(n,1))

    # test set inclusion
    for k in range(n):
        v = p[k]  # next vertex to be eliminated

        # indices of neighbors that correspond to strictly lower triangular elements in reordered pattern
        r = set([rj for rj in ri[cp[v]:cp[v+1]] if ip[rj] > k])  

        for rj in r:
            if not r.issubset(set(ri[cp[rj]:cp[rj+1]])): return False
            
    return True

def merge_size_fill(tsize = 8, tfill = 8):
    """
    Simple heuristic for supernodal amalgamation (clique
    merging). 

    Returns a function that returns `True` if either (i) supernode k and
    supernode par(k) are both of size at most `tsize`, or (ii),
    merging supernodes par[k] and k induces at most `tfill` nonzeros
    in the lower triangular part of the sparsity pattern.

    :param tsize:   nonnegative integer; threshold for merging based on supernode sizes
    :param tfill:   nonnegative integer; threshold for merging based on induced fill
    """
    assert tsize >= 0, "tsize must be nonnegative"
    assert tfill >= 0, "tfill must be nonnegative"
    
    def fmerge(ccp,cck,np,nk):
        """
        Supernodal merge heuristic.

           d = fmerge(Jp, Jk, Np, Nk)

        PURPOSE
        Returns true if either (i) supernode k and supernode par(k)
        are both of size at most %i, or (ii), merging supernodes
        par(k) and k induces at most %i edges.

        ARGUMENTS
        Jp        integer; size of parent clique

        Jk        integer; size of child clique

        Np        integer; size of parent supernode 

        Nk        integer; size of child supernode
         
        RETURNS
        d         boolean
        """
        fill = (ccp - (cck - nk))*nk
        if fill <= tfill or (nk <= tsize and np <= tsize):
            return True
        else:
            return False
    # insert parameters into docstring of fmerge and return fmerge
    fmerge.__doc__ %= (tsize, tfill)
    return fmerge

class symbolic(object):
    """
    Symbolic factorization object.

    Computes symbolic factorization of a square sparse matrix
    :math:`A` and creates a symbolic factorization object.

    :param A:   :py:class:`spmatrix`
    :param p:   permutation vector or ordering routine (optional) 
    :param merge_function:  routine that implements a merge heuristic (optional)

    The optional argument `p` can be either a permutation vector or an
    ordering rutine that takes an :py:class:`spmatrix` and returns a
    permutation vector.

    The optional argument `merge_function` allows the user to
    merge supernodes in the elimination tree in a greedy manner;
    the argument must be a routine that takes the following four arguments and
    returns either `True` or `False`:
	
    :param cp:   clique order of the parent of clique :math:`k`
    :param ck:   clique order of clique :math:`k`
    :param np:   supernode order of the parent of supernode :math:`k`
    :param nk:   supernode order of supernode :math:`k`

    The clique k is merged with its parent if the return value is `True`.
    """
    
    def __init__(self, A, p = None, merge_function = None):

        assert isinstance(A,spmatrix), "A must be a sparse matrix"
        assert A.size[0] == A.size[1], "A must be a square matrix"

        # Symmetrize A
        Ap = symmetrize(A)
        nnz_Ap = (len(Ap)+Ap.size[0])/2   # number of lower triangular nonzeros

        # Permute if permutation vector p or ordering routine is specified
        if p is not None:
            if type(p) is BuiltinFunctionType or type(p) is FunctionType:
                p = p(Ap)
            assert len(p) == A.size[0], "length of permutation vector must be equal to the order of A"
            Ap = perm(Ap,p)
            ip = matrix(0,(A.size[0],1))
            ip[p] = matrix(range(len(ip)))
        else:
            ip = p
            
        self.__p = p
        self.__ip = ip

        # Symbolic factorization
        par = etree(Ap)
        post = post_order(par)
        colcount = counts(Ap, par, post)
        nnz_Ae = sum(colcount)
        snode, snptr, snpar = supernodes(par, post, colcount)
        snpost = post_order(snpar)
        if merge_function:
            colcount, snode, snptr, snpar, snpost = amalgamate(colcount, snode, snptr, snpar, snpost, merge_function)
        sncolptr, snrowidx = embed(Ap, colcount, snode, snptr, snpar, snpost)
        relptr, relidx = relative_idx(sncolptr, snrowidx, snptr, snpar)
        
        # build chptr
        chptr = matrix(0, (len(snpar)+1,1))
        for j in snpost: 
            if snpar[j] != -1: chptr[snpar[j]+1] += 1
        for j in range(1,len(chptr)):
            chptr[j] += chptr[j-1]

        # build chidx
        tmp = +chptr
        chidx = matrix(0,(chptr[-1],1))
        for j in snpost:
            if snpar[j] != -1: 
                chidx[tmp[snpar[j]]] = j
                tmp[snpar[j]] += 1
        del tmp

        # compute stack size
        stack_size = 0
        stack_mem = 0
        stack_max = 0
        frontal_max = 0
        stack = []
        for k in snpost:
            nn = snptr[k+1]-snptr[k]
            na = relptr[k+1]-relptr[k]
            nj = na + nn
            frontal_max = max(frontal_max, nj**2)
            for j in range(chptr[k+1]-1,chptr[k]-1,-1):
                v = stack.pop()
                stack_mem -= v**2
            if (na > 0):
                stack.append(na)
                stack_mem += na**2
                stack_max = max(stack_max,stack_mem)
                stack_size = max(stack_size,len(stack))
        self.frontal_len = frontal_max
        self.stack_len = stack_max
        self.stack_size = stack_size
                
        # build blkptr
        blkptr = matrix(0, (len(snpar)+1,1))
        for i in range(len(snpar)):
            blkptr[i+1] = blkptr[i] + (snptr[i+1]-snptr[i])*(sncolptr[i+1]-sncolptr[i])

        self.__n = len(snode)
        self.__Nsn = len(snpar)
        self.__snode = snode
        self.__snptr = snptr
        self.__chptr = chptr
        self.__chidx = chidx
        self.__snpar = snpar
        self.__snpost = snpost
        self.__relptr = relptr
        self.__relidx = relidx
        self.__sncolptr = sncolptr
        self.__snrowidx = snrowidx
        self.__blkptr = blkptr
        self.__fill = (nnz_Ae-nnz_Ap,self.nnz-nnz_Ae)

        return

    def __repr__(self):
        return "<symbolic factorization, n=%i, nnz=%i, nsn=%i>" % (self.n,self.nnz,self.Nsn)

    def __str__(self):
        opts = printing.options
        printing.options = {'iformat':'%1i','dformat':'%1.0f',\
                            'width':printing.options['width'],'height':printing.options['height']}
        tmp = self.sparsity_pattern(reordered = True, symmetric = True).__str__()
        printing.options = opts
        return tmp.replace('0',' ').replace('1','X')

    def sparsity_pattern(self, reordered = True, symmetric = True):
        """
        Returns a sparse matrix with the filled pattern. By default,
        the routine uses the reordered pattern, and the inverse
        permutation is applied if `reordered` is `False`.
        
        :param reordered:  boolean (default: `True`)
        :param symmetric:  boolean (default: `True`)	
        """
        return cspmatrix(self, 1.0).spmatrix(reordered = reordered, symmetric = symmetric)

    @property
    def n(self):
        """Number of nodes (matrix order)"""
        return self.__n

    @property
    def Nsn(self):
        """Number of supernodes"""
        return self.__Nsn

    @property
    def snode(self):
        """
        Supernode array: supernode :math:`k` consists of nodes
        `snode[snptr[k]:snptr[k+1]]` where `snptr` is the supernode
        pointer array
        """
        return self.__snode

    @property
    def snptr(self):
        """
        Supernode pointer array: supernode :math:`k` is of order
        `snpptr[k+1]-snptr[k]` and supernode :math:`k` consists of nodes
        `snode[snptr[k]:snptr[k+1]]`
        """
        return self.__snptr

    @property
    def chptr(self):
        """
        Pointer array associated with `chidx`:
        `chidx[chptr[k]:chptr[k+1]]` are the indices of the children
        of supernode k.
        """
        return self.__chptr

    @property
    def chidx(self):
        """
        Integer array with indices of child vertices in etree: 
        `chidx[chptr[k]:chptr[k+1]]` are the indices of the children
        of supernode :math:`k`.
        """
        return self.__chidx


    @property
    def snpar(self):
        """
        Supernode parent array: supernode :math:`k` is a root of the
        supernodal elimination tree if `snpar[k]` is equal to -1, and
        otherwise `snpar[k]` is the index of the parent of supernode
        :math:`k` in the supernodal elimination tree
        """
        return self.__snpar

    @property
    def snpost(self):
        """Supernode post-ordering"""
        return self.__snpost

    @property
    def relptr(self):
        """Pointer array assoicated with `relidx`."""
        return self.__relptr

    @property 
    def relidx(self): 
        """ The relative index array facilitates fast "extend-add" and
        "extract" operations in the supernodal-multifrontal
        algorithms. The relative indices associated with supernode
        :math:`k` is a list of indices :math:`I` such that the frontal
        matrix :math:`F` associated with the parent of node :math:`k`
        can be updated as `F[I,I] += Uj`. The relative indices are
        stored in an integer array `relidx` with an associated pointer
        array `relptr`."""
        return self.__relidx

    @property
    def sncolptr(self):
        """
        Pointer array associated with `snrowidx`.
        """
        return self.__sncolptr

    @property
    def snrowidx(self):
        """
        Row indices associated with representative vertices:
        `snrowidx[sncolptr[k]:sncolptr[k+1]]` are the row indices in
        the column corresponding the the representative vertex of
        supernode :math:`k`, or equivalently,
        `snrowidx[sncolptr[k]:sncolptr[k+1]]` is the :math:`k`'th
        clique.
        """
        return self.__snrowidx

    @property
    def blkptr(self):
        """
        Pointer array for block storage of chordal sparse matrix.
        """        
        return self.__blkptr

    @property
    def fill(self):
        """
        Tuple with number of lower-triangular fill edges: `fill[0]` is
        the fill due to symbolic factorization, and `fill[1]` is the
        fill due to supernodal amalgamation"""
        return self.__fill
    
    @property
    def nnz(self):
        """
        Returns the number of lower-triangular nonzeros.
        """
        nnz = 0
        for k in range(len(self.snpost)):
            nn = self.snptr[k+1]-self.snptr[k]    
            na = self.relptr[k+1]-self.relptr[k]
            nnz += nn*(nn+1)/2 + nn*na
        return nnz

    @property
    def p(self):
        """
        Permutation vector
        """
        return self.__p

    @property
    def ip(self):
        """
        Inverse permutation vector
        """
        return self.__ip

    @property
    def clique_number(self):
        """
        The clique number (the order of the largest clique)
        """
        return max([self.sncolptr[k+1]-self.sncolptr[k] for k in range(self.Nsn)])

    
    def cliques(self):
        """
        Returns a list of cliques (reordered pattern)
        """
        return [list(self.snrowidx[self.sncolptr[k]:self.sncolptr[k+1]]) for k in range(self.Nsn)]

    def separators(self):
        """
        Returns a list of separator sets (reordered pattern)
        """
        return [list(self.snrowidx[self.sncolptr[k]+self.snptr[k+1]-self.snptr[k]:self.sncolptr[k+1]]) for k in range(self.Nsn)] 

    def supernodes(self):
        """
        Returns a list of supernode sets (reordered pattern)
        """
        return [list(self.snode[self.snptr[k]:self.snptr[k+1]]) for k in range(self.Nsn)]

    def parent(self):
        """
        Returns a supernodal parent list: the i'th element is equal to -1 if 
        supernode i is a root node in the clique forest, and otherwise
        the i'th element is the index of the parent of supernode i.
        """
        return list(self.snpar)
    

