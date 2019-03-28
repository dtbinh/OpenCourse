% Main script for testing various methods.
% The algorithms in these files are explained in more detail in:
% http://www.seas.ucla.edu/~vandenbe/publications/spcones.pdf
% 
%
% For a sparsity pattern E, the three problems are
%
%   1) the sparse SDP cone (sdp)
%
%       minimize    ||X-A||_F^2
%       subject to  X >= 0,
%                   X in sparsity pattern E
%
%   2) the cone of sparse matrices with SDP completion (sdpc)
%
%       minimize ||X-A||_F^2
%       subject to X in sparsity pattern E,
%                  X has a PSD completion
%
%   3) the cone of sparse matrices with EDM completion (edm)
%
%       minimize ||X-A||_F^2
%       subject to X in sparsity pattern E,
%                  X has an EDM completion
%
%
% To use an interior point method, assign "method" to the appropriate cvx
% solver (e.g. Sedumi, sdpt3, mosek). Both CVX and this solver must 
% already be installed.
%
% The available first-order chordal decomposition solvers are 
%   1) the proximal gradient method and the fast proximal gradient method
%     (chordal sparsity pattern and SDP cones only)
%   1) block coordinate ascent
%     (chordal sparsity pattern and SDP cones only)
%   2) the Douglas-Rachford method
%
% Authors: Yifan Sun & Lieven Vandenberghe 
% Date: March 2015

clc
clear
addpath(genpath('.'))

% If IPMs desired, assign correct path to CVX. Otherwise, comment out this
% block.
addpath(genpath(path_to_cvx))
cvx_setup



%% Parameters

%%% ENTER CHOICE OF METHODS AND PROBLEMS HERE %%%
% String label for methods:
%   'sedumi', 'sdpt3', 'mosek': corresponding IPM
%   'PG', 'BCD', 'PDR', 'DDR' : decomposition method
% enter the solvers for comparison in list method_list
%
% String label for problems:
% 'sdp'  : sparse positive semidefinite cone
% 'sdpc' : sparse positive semidefinite completable cone
% 'edmc'  : sparse Euclidean distance matrix completable cone
%
% Boolean indicators
%   ischordal : 1 if cliques are for chordal sparsity pattern 
%               0 if cliques are for chordal embedding of sparsity pattern



%Example

problemtype = 'edmc';
method_list = [{'sdpt3'},{'PG'},{'BCD'},{'PDR'},{'DDR'}];
ischordal = true;

%% Display info
if ischordal
    fprintf('Projection onto sparse (chordal) %s cone\nMethods:  ', problemtype )
else
    fprintf('Projection onto sparse (nonchordal) %s cone\nMethods: ', problemtype)
end
for ml = 1:length(method_list)-1
    fprintf('%s, ',method_list{ml})
end
fprintf('and %s.',method_list{end})

%% Generate problem


%%%CHOOSE ONE OF THE NEXT TWO METHODS TO GENERATE SPARSITY PATTERN
%load a matrix
load('can_61.mat')
N = size(Problem.A,1);
Sporig = sparse(Problem.A);
Sporig = Sporig + Sporig' + speye(N);
Sporig = Sporig > 0;
A = Sporig  .* randn(N);
A = (A + A')/2;
A = sparse(A);


% sensor network localization
% N = 25;
% u = rand(N,3);
% A = u*u';
% A = diag(A)*ones(1,N) + ones(N,1)*diag(A)' - 2*A;
% rho = 10; 
% Sporig = A < (rho/N);
% Sporig = Sporig + Sporig' + speye(N);
% Sporig = Sporig > 0;
% noise = randn(N) * 0.01;
% noise = (noise + noise') / 2;
% A = (A.*Sporig) + noise.*Sporig;

%% Chordal embedding
%%% MODIFY PROBLEM PARAMETERS HERE %%%
tsize = 5; tfill = 5; 
% Sp is permuted original sparsity, Spfill is permuted and filled
[Sp, Spfill,Sp_unmerged, p, overlap, veclen, cliques, cliques_unmerged] = chordal_embedding(Sporig,tfill,tsize, true); 
L = max_overlap(cliques);


% Additional preprocessing
Ap = A(p,p);
if ischordal
    D = Spfill;
else
    D = Sp;
end


figure(1)
clf
subplot(2,2,1)
spy(Sporig)
title('Original sparsity pattern')
subplot(2,2,2)
spy(Sp)
title('Permuted sparsity pattern')
subplot(2,2,3)
spy(Sp_unmerged)
title('Filled sparsity pattern')
subplot(2,2,4)
spy(Spfill)
title('Merged sparsity pattern')


%% Interior point methods

%Using Sedumi solver
if sum(strcmpi(method_list,'sedumi'))
    [Xipm, out] = IPM('sedumi',problemtype, Ap, D);
    fprintf('Obj: %f, CPU time = %f\n\n',out.obj, out.runtime)
end

%Using SDPT3 solver
if sum(strcmpi(method_list,'sdpt3'))
    [Xipm, out] = IPM('sdpt3',problemtype, Ap, D);
    fprintf('Obj: %f, CPU time = %f\n\n',out.obj, out.runtime)
end

%Using Mosec solver
if sum(strcmpi(method_list,'mosek'))
    [Xipm, out] = IPM('mosek',problemtype, Ap, D);
    fprintf('Obj: %f, CPU time = %f\n\n',out.obj, out.runtime)
end


%% options for first-order algorithms.
opts.cliques = cliques;             % list of cliques stored as cells
opts.problemtype = problemtype;     % sdp for P1, sdpc for P2, edmc for P3
opts.D = D;                         % sparsity pattern E of problem
opts.verbose = true;                % allow printout
opts.epoch = 25;                    % # iterations between checking stopping conditions
opts.maxiter = 1000;                % total allowed number of iterations
opts.tol = 1e-3;                    % tolerance for stopping criteria
opts.maxtime = 60*60*4;             % total allowed number of seconds

opts = preprocess_opts(A,opts);     % automatically generate remaining parameters 


%% Chordal decomposition methods

%Fast proximal gradient (for chordal problems only)
if sum(strcmpi(method_list,'PG'))
    if ~ischordal
        fprintf('Cannot run prox gradient for nonchordal E\n')
    else
        opts.accelerate = true;
        opts.L = L;
        [Xfo,out] = proxgrad(Ap, opts);
        fprintf('Obj: %f, CPU time = %f, res: %f\n\n',out.obj(end), out.runtime(end), out.err(end).primal)
    end
end

%Block coordinate descent / Dykstra's method (for chordal problems only)
if sum(strcmpi(method_list,'BCD'))
    if ~ischordal
        fprintf('Cannot run block coordinate descent for nonchordal E\n')
    else
        [Xfo, out] = bcd(Ap, opts);
        fprintf('Obj: %f, CPU time = %f, res: %f\n\n',out.obj(end), out.runtime(end), out.err(end).primal)
    end
end

%Douglas-Rachford primal problem formulation
if sum(strcmpi(method_list,'PDR'))
    opts.t = 1;
    opts.rho = 1.75;
    opts.isdual = false;
    [Xfo, out] = dr(Ap, opts);
    fprintf('Obj: %f, CPU time = %f, res: (%f,%f) \n\n',out.obj(end), out.runtime(end),out.err(end).primal,out.err(end).dual)
end
            
%Douglas-Rachford dual problem formulation
if sum(strcmpi(method_list,'DDR'))
    opts.t = 1;
    opts.rho = 1.75;
    opts.isdual = true;
    [Xfo, out] = dr(Ap, opts);
    fprintf('Obj: %f, CPU time = %f, res: (%f,%f) \n\n',out.obj(end), out.runtime(end),out.err(end).primal,out.err(end).dual)
end



        
