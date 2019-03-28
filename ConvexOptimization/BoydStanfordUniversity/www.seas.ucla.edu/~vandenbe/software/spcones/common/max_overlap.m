function L = max_overlap(cliques)
    % MAX_OVERLAP : Find maximum number of times an index is overlapped in a set of
    % cliques
    %
    %   L = max_overlap(cliques) returns the maximum overlap amongst cliques. 
    %
    % INPUT: 
    %   cliques : a cellarray containing the indices beta_k of the cliques
    %
    % OUTPUT:
    %   L : max overlap constant
    %
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015

    % The primary use of this function is to determine the Lipschitz
    % constant in the proximal gradient methods, for the smooth objective
    % term :
    %     f(Z) = norm(sum_k Pk.T * Zk* Pk - A ,'fro')^2
    %
   
    I = [];
    for k = 1:length(cliques)
        I = [I ; cliques{k}];
    end
    Lvec = sparse(I,1,1);
    L = max(Lvec);
end
