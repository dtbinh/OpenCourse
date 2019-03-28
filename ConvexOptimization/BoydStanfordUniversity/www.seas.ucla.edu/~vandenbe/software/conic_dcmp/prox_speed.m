% Collects runtime of prox evaluations using CVX and SDPT3 for dense cliques.
% Solves
% 
%     minimize trace(C*X) + 1/2 * ||X-Z||_F^2
%     subject to
%         A(X) == b
%         X >= 0
%         
%         
% Requires problems to be first generated (via generate_problem.py) to maintain 
% consistency with fast prox results.

clc
clear
cvx_setup

vec = 25:25:250;

for pin = 1:length(vec)
    
    p = vec(pin); % order of matrix
    s = p; % number of constraints
    
    % load precomputed data
    fi = sprintf('problems/p%d_s%d.mat',p,s)
    try
        load(fi,'C', 'A', 'b', 'Z')
        C; A; b; Z;
    catch
        error(sprintf('Problem p = %d, s = %d not yet generated. Please run generate_problem.py',p,s))
    end
    
    for solverid = 1:4
        fprintf('p = %d, solver %d\n',p,solverid)

        if solverid == 1
            %cvx, sdpt3
            [t,o,i] = solve_single_prox_via_cvx(A,b,C,Z, 'sdpt3');
         
        elseif solverid == 2
            %direct sdpt3 using qcone, one big socp
            [t,o,i,AA,bb,cc,K] = solve_single_prox_via_ipm(A,b,C,Z, 1, 'q',0);
            [blk,At,Ct,bt] = read_sedumi(AA,bb,cc,K);
            [o,X,y,S, info] = sqlp(blk,At,Ct,bt,sqlparameters);
            t = info.cputime;
            o = o(1);
            i = info.iter;
            
        elseif solverid == 3
            %direct sdpt3 using qcone, p little SOCPs for each column
            [t,o,i,AA,bb,cc,K] = solve_single_prox_via_ipm(A,b,C,Z, p, 'q',0);
            [blk,At,Ct,bt] = read_sedumi(AA,bb,cc,K);
            
            [o,X,y,S, info] = sqlp(blk,At,Ct,bt,sqlparameters);
            t = info.cputime;
            o = o(1);
            i = info.iter;
        elseif solverid == 4
            %direct sdpt3 using qcone, p^2 little SOCPs for each i,j
            [t,o,i,AA,bb,cc,K] = solve_single_prox_via_ipm(A,b,C,Z, p^2, 'q',0);
            [blk,At,Ct,bt] = read_sedumi(AA,bb,cc,K);
            [o,X,y,S, info] = sqlp(blk,At,Ct,bt,sqlparameters);
            t = info.cputime;
            o = o(1);
            i = info.iter;
            
        end
        
        obj(solverid,pin) = o;
        time(solverid,pin) = t;
        itr(solverid,pin) = i;
        
    end
    
end

figure(1)
clf
plot(time)