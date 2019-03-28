from cvxopt import matrix,spmatrix, sparse,amd,spdiag
from symbolic import symbolic, merge_size_fill
import pylab
import numpy as np
from matplotlib.pyplot import *
import matfile
import pickle
from misc import reshape, spdiag_rect


        
def convert(Aun, bun, Cun, S, dims, Atype = spmatrix, tfill = 0, tsize = 0):
    """
    Implements conversion where block-diagonal correlative sparsity is possible.
    
    Given an SDP problem with symmetric matrix variable X
    
    minimize trace(Cun*X)
    subject to
        Aun.T*X[:] = bun
        X == positive semidefinite
        
    converts to the problem with matrix variables Xi
    
    minimize sum_i trace(Ccon[i]*Xi)
    subject to
        Acon[i].T*Xi = bcon[i] for all i
        Xi == positive semidefinite for all i
    
    ARGUMENTS:
        Aun, bun, Cun : coefficients of unconverted problem
        S : sparsity pattern of unconverted problem 
            (Sij = 1 if (i,j) in aggregate sparsity pattern, Sij = 0 otherwise)
        tfill, tsize : parameters in Acon = Acon['s']
Ccon = Ccon['s']clique merging
        
    RETURNS:
        Acon, bcon, Ccon : list of coefficient matricies for each clique
        Lcon : sparse matrix containing consistency constraints (L.T*Xt = 0 
                where Xt is a stacked vector containing Xi's)
        dims: dictionary of dimensions:
            dims['l'] : number of linear variables in nonnegative cone
            dims['q'] : list of sizes of quadratic cones
            dims['s'] : list of orders of semidefinite cones
    """
    n = dims['s'][0]
    nl = dims['l']
    m = Aun.size[1]

    #Find clique tree decomposition of sparsity pattern
    if tsize > 0 or tfill > 0:
        symb = symbolic(S,  p = amd.order, merge_function = merge_size_fill(tsize = tsize, tfill = tfill))
    else:
        symb = symbolic(S,  p = amd.order)
    ns = [len(clique) for clique in symb.cliques()]
    
    #permute coefficient and sparsity matricies
    Cunl = Cun[:nl]
    Cuns = reshape(Cun[nl:],n,n)
    Cuns = +Cuns[symb.p,symb.p]
    S = S[symb.p,symb.p]
    cliques = symb.cliques()
    
    
    Aunl = Aun[:nl,:]
    Auns = Aun[nl:,:]
    I,J,v = Auns.I,Auns.J,Auns.V
    j = [int(k/n) for k in list(I)]
    i = [int(I[k]-(j[k]*n)) for k in xrange(len(I))]
    Auns = spmatrix(v,symb.ip[i] + symb.ip[j]*n,J,Aun.size)

    cliques_all = []
    for cl in cliques: cliques_all.extend(cl)
    assert n == len(set(cliques_all)), "Each index must be included in at least one clique."

    #Convert C
    Ccon = []
    for clique in cliques:
        Ccon.append(+Cuns[clique,clique])
        Cuns[clique,clique] = 0.0


    #Convert A. Greedily assign constraints to cliques in which 
    #nonzeros of Ak are fully contained in clique k
    
    Aconl = [[matrix(0.0,(0,nl))] for i in ns]
    Acons = [[matrix(0.0,(0,i**2))] for i in ns]
    bcon = [[] for k in xrange(len(ns))]
    
    constraint_order =[] 
    for k in xrange(len(ns)):
        clique = cliques[k]
        unwrapped = [(i+j*n) for i in list(clique) for j in list(clique)]
        for con in xrange(m):
            if con in constraint_order: continue
            AI = sparse(Auns[:,con]).I
            if set(AI).issubset(set(unwrapped)):
                Aconl[k].append(Aunl[:,con].T)
                Acons[k].append(Auns[unwrapped,con].T)
                bcon[k].append(bun[con])
                constraint_order.append(con)
        if Atype is matrix: 
            Aconl[k] = matrix(Aconl[k]).T  
            Acons[k] = matrix(Acons[k]).T   
        else: 
            Aconl[k] = sparse(Aconl[k]).T  
            Acons[k] = sparse(Acons[k]).T
        
        bcon[k] = matrix(bcon[k])
    assert len(constraint_order) == Aun.size[1], "Conversion for not totally seperable not yet implemented"
    
    Acon = sparse([sparse([a.T for a in Aconl]).T,spdiag_rect(Acons)])
    bcon = matrix(bcon)
    Ccon = matrix([Cunl,matrix([c[:] for c in Ccon])])

    
    #get consistency constraints L*x = 0
    total_variables = [0]
    for k in ns:
        total_variables.append(total_variables[-1] + k**2)
        
    L = ([],[],[])
    ncol = 0
    for child in xrange(len(symb.snpar)):
        parent = symb.snpar[child]
        if parent == -1: continue

        for ip in xrange(len(cliques[parent])):
            for jp in xrange(ip+1):
              
                for ic in xrange(len(cliques[child])):
                    if cliques[child][ic] != cliques[parent][ip]: continue
                    for jc in xrange(ic+1):
                        if cliques[child][jc] != cliques[parent][jp]: continue
                        
                        
                        zp = jp * len(cliques[parent]) + ip + total_variables[parent]
                        zc = jc * len(cliques[child]) + ic + total_variables[child]


                        L[0].extend([1.0,-1.0])
                        L[1].extend([zp,zc])
                        L[2].extend([ncol,ncol])
                        
                        
                        if jp != ip:
                            zp = ip * len(cliques[parent]) + jp + total_variables[parent]
                            zc = ic * len(cliques[child]) + jc + total_variables[child]
                        
                            L[0].extend([1.0,-1.0])
                            L[1].extend([zp,zc])
                            L[2].extend([ncol,ncol])
                        ncol += 1

    Lcon = spmatrix(L[0],matrix(L[1])+nl,L[2],size=(total_variables[-1]+nl,ncol))

    dims = {'l':nl,'q':[],'s':ns}
    return Acon,bcon,Ccon,Lcon, dims
