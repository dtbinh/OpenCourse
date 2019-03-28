function x = project_over_cliques(project, x,cliques)
    % PROJECT_OVER_CLIQUES : Performs a lower dimensional projection on 
    % each clique of a sparse matrix in vectorized form.
    %
    %     x = project_over_cliques(project, x,cliques): Given a function 
    %     pointer project(x) that performs the desired projection, this 
    %     function projects on each vectorized submatrix X_i
    %
    % INPUTS
    %   project : function pointer for projection
    %   x       : stacked vector (x_1,...,x_l) for lower dimensional
    %             vectors x_k
    %   cliques : a cellarray containing the indices beta_k of the cliques
    %
    % OUTPUTS
    %   x       : stacked vector (xt_1,...,xt_l) where 
    %             x_k = project(x_k)
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015

    offset = 0;
    for k = 1:length(cliques)      
        cl = cliques{k};
        ncl = length(cl);
        tmp = project(reshape(x(offset + 1:offset + ncl^2),ncl,ncl));
        tmp = vec(tmp);
        x(offset + 1:offset + ncl^2) = tmp;
        offset = offset + ncl^2;
        
    end
end
