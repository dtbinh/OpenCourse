function [time,obj,itr,Ased,bsed,csed,K,X,Y] = solve_single_prox_via_ipm(A,b,C,Z, ncones, cone,tosolve)
    
    % Solves
    % 
    %     minimize trace(C*X) + 1/2 * ||X-Z||_F^2
    %     subject to
    %         A(X) == b
    %         X >= 0
    %
    % using SDPT3. The term ||X-Z||_F^2 can be modeled using quadratic
    % cones in a number of ways. 

    time = -1; obj = -1; itr = -1;
    p = size(C,1);
    
    nlower = p*(p+1)/2;
    noff = p*(p-1)/2;
    tmp = eye(p);
    diag_index = find(tmp(:));
    tmp = tril(ones(p))-eye(p);
    offlower_index = find(tmp(:));
    M = reshape(1:p^2,p,p);
    M = M - tril(M);
    tmp = M';
    offupper_index= tmp(:);
    offupper_index = offupper_index(offupper_index>0);

    if (ncones == 1) && strcmp(cone, 'q')
        
        
        % minimize trace(C*X) + 1/2*||X-Z||^2
        
        K.f = 1;
        K.q = nlower+2;
        K.s = p;
        ms = size(A,2);
        csed = [1/2 ;              zeros(nlower+2,1) ;       C(:) ];
        Ased = [[spalloc(ms,1,0) ,   spalloc(ms,nlower+2,0) ,	A' ];...
            [-1 ,                  1, spalloc(1,nlower+1,0) ,    spalloc(1,p^2,0) ]; ...
            [-1,                   spalloc(1,nlower+1,0), 1,     spalloc(1,p^2,0)];...
            [spalloc(p,2,0),         speye(p),spalloc(p,noff + 1,0),  -2*sparse(1:p,diag_index,1,p,p^2,p)]; ...
            [spalloc(noff,2,0),      spalloc(noff,p,0), speye(noff), spalloc(noff,1,0),    -2*sqrt(2)/2*sparse([1:noff,1:noff],[offlower_index,offupper_index],1,noff, p^2, 2*noff)]];
        bsed = [b ; 1 ; -1 ; -2*diag(Z) ; -2*sqrt(2)*Z(offlower_index)];
        if tosolve
            [X,Y,INFO] = sedumi(Ased,bsed,csed,K);
            time = INFO.cpusec;
            obj = csed'*X;
            itr = INFO.iter;
            
        end
        
        
        
    elseif (ncones == p^2) && strcmp(cone, 'q')
        
        
        % minimize trace(C*X) + 1/2*sum_{ij} (Xij-Zij)^2
        
        
        
        
        K.f = nlower;
        K.q = ones(nlower,1)*3;
        K.s = p;
        ms = size(A,2);
        csed = [1/2*ones(nlower,1) ;  zeros(3*nlower,1) ; C(:) ];
        Ased = [spalloc(ms,nlower,0) ,  spalloc(ms,3*nlower,0),          A' ];
        
        bsed = b;
        
        [I,J,V] = find(Ased);
     
        Ioff = ms;
        for j = 1:p
            
            
            k = j+(j-1)*p;
            index = p*(j-1)+j-(j-1)*j/2;
            I = [I ; [1,3, 1,2,3, 2,]'+Ioff];
            J = [J ; [index, index, nlower+(index-1)*3+1, nlower+(index-1)*3+2,nlower+(index-1)*3+3,nlower+nlower*3+k]'];
            V = [V ; [-1,-1, 1,1,1,-2]'];

            Ioff = Ioff + 3;
            bsed = [bsed ; 1; -2*Z(j,j); -1];

            if j < p
                k = j+(j-1)*p+(1:p-j);
                k2 = j+(j+(1:p-j)-1)*p;
                off = index + (1:(p-j));
                
                tmp =  repmat([1,3, 1,2,3, 2,2]',1,p-j) + repmat(0:3:3*(p-j-1),7,1);
                I = [I;Ioff + tmp(:)];
               Ioff = Ioff + 3*(p-j);
               tmp = [off; off; nlower+(off-1)*3+1; ...
                    nlower+(off-1)*3+2;nlower+(off-1)*3+3;...
                    nlower+nlower*3+k;...
                    nlower+nlower*3+k2];
                J =[J;tmp(:)];
                V = [V;repmat([-1,-1, 1,1,1,-2/2*sqrt(2),-2/2*sqrt(2)],1,p-j)'];

                tmp = [ ones(1,p-j) ; -2*sqrt(2)*Z((j+1):end,j)' ; -ones(1,p-j)];
                bsed = [bsed ; tmp(:)];
            end
        end
        Ased = sparse(I,J,V);
        
        if tosolve
            [X,Y,INFO] = sedumi(Ased,bsed,csed,K);
            time = INFO.cpusec;
            obj = csed'*X;
            itr = INFO.iter;
            
        end
   
        
    elseif (ncones == p) && strcmp(cone, 'q')
        
        

        % minimize trace(C*X) + 1/2*sum_{i} ||X(i,:)-Z(i,:)||^2
        
        K.f = p;
        K.q = (p:-1:1)+2;
        K.s = p;
        
        ms = size(A,2);
        csed = [1/2*ones(p,1) ;   zeros(nlower + 2*p,1) ;     C(:) ];
        Ased = [spalloc(ms,p,0) , spalloc(ms,nlower + 2*p,0),   A' ];
        bsed = b;
        
        [I,J,V] = find(Ased);
        Ioff = ms;
        offset = 0;
        
        for j = 1:p
            cl = (p-j);

            tmp = [[sparse(1,j,-1,1,p,1), sparse(1,offset + 1,1,1,2*p+nlower,1),  spalloc(1,p^2,0)];...
            [sparse(1,j,-1,1,p,1), sparse(1,offset + 2,1,1,2*p+nlower,1),  spalloc(1,p^2,0)]];
            [ii,jj,vv] = find(tmp);
           
            I = [I ; ii+Ioff]; J = [J ; jj]; V = [V ; vv];
             Ioff = Ioff + 2;
            
            bsed = [bsed ; 1; -1];

            bsed = [bsed ; -2*Z(j,j)];
            
            tmp = [spalloc(1,p,0), sparse(1,offset+2+1,1,1,2*p+nlower, 1),   sparse(1,(j-1)*p+j,-2,1,p^2, 1)];
            [ii,jj,vv] = find(tmp);
            
            I = [I ; ii'+Ioff]; J = [J ; jj']; V = [V ; vv'];
            Ioff = Ioff + 1;
           
            if j < p
                tmp = [spalloc(cl,p,0), sparse(1:cl,offset+2+(2:(cl+1)),1,cl,2*p+nlower, cl),   sparse([1:cl,1:cl],[(j-1)*p+j + (1:cl)*p,(j-1)*p+j+(1:cl)],-2/2*sqrt(2),cl,p^2, 2*cl)];
                [ii,jj,vv] = find(tmp);
             
                
                if size(ii,2) > size(ii,1)
                    ii = ii'; jj = jj'; vv = vv';
                end
                I = [I ; ii+Ioff]; J = [J ; jj]; V = [V ; vv];
                Ioff = Ioff + cl;
                bsed = [bsed ; -2*sqrt(2)*Z((j+1):end,j)];
                
            end
            offset = offset + 2+(p-j+1);
        end
        Ased = sparse(I,J,V);
        if tosolve
            [X,Y,INFO] = sedumi(Ased,bsed,csed,K);
            time = INFO.cpusec;
            obj = csed'*X;
            itr = INFO.iter;
        end
        
    end
end
