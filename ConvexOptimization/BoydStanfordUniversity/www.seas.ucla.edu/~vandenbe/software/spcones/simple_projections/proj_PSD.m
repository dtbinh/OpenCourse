function X = proj_PSD(X,sp)
    % PROJ_PSD : Projection onto the dense PSD cone
    %
    %   X = proj_PSD(X, sp) : projects matrix X onto dense positive
    %   semidefinite cone.
    %
    % INPUT: 
    %   X   : N x N matrix
    %   sp  : In the eigenvalue decomposition, 1 to use eigs (Lanczos-type), 
    %         0 to use eig (full eigenvalue decomposition)
    % OUTPUT:
    %   X   : projected N x N matrix
    %
    % Authors: Yifan Sun & Lieven Vandenberghe 
    % Date: March 2015

    N = size(X,1);
    if N == 1
        X = max(X, 0);
    end
    X = .5*(X+X');
    if sp
        
    	[~,d,~] = ldl(X);
        num_neig = full(sum(diag(d)<-0));
        num_peig = full(sum(diag(d)>0));

        if num_neig < num_peig
            [Vvec,lambda] = eigs(X,num_neig,'sa');
            X = X - Vvec*lambda*Vvec';
        else
            [Vvec,lambda] = eigs(X,num_peig,'la');
            X =  Vvec*lambda*Vvec';
        end
    else
        if issparse(X)
            X = full(X);
        end
        
        flag = true;
        addon = 0;
        while(flag)
            try
                [Vvec,lambda] = eig(X + addon * eye(size(X,1)));
                flag = false;
            catch
                addon = addon + 1;
            end
        end
        lambda = diag(lambda) - addon;
        lambda(lambda<0) = 0;
        
        
        if length(lambda) > 5000
            VL = Vvec * 0;
            block = 1000;
            offset = 0;
            while offset < length(lambda)
                block = min(block, length(lambda) - offset);
                index = offset + (1:block);
                VL(:,index) = VL(:,index)*diag(lambda(index));
                offset = offset + block;
            end
        else
            VL = Vvec * diag(lambda);
        end
        X = VL*Vvec';
    end
 
    X = .5*(X+X');
end
