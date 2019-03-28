function [x,status,lambda,nu] = gpposy(A,b,szs,varargin)

% --------------------------------------------------------------------------------
% gpposy: A Matlab solver for Geometric Programs (GP) in posynomial form
%
% Kwangmoo Koh, Seung-Jean Kim, Almir Mutapcic, and Stephen Boyd. (2006)
% gpposy is part of the ggplab package: http://www.stanford.edu/~boyd/ggplab/
% --------------------------------------------------------------------------------
%
% [x,status,lambda,nu] = gpposy(A,b,szs,G,h,l,u,quiet)
%
% solves the geometric program in posynomial form
%
%  minimize    sum_k b0_k*x1^A0_{k1}*...*xn^A0_{kn}
%  subject to  sum_k bi_k*x1^Ai_{k1}*...*xn^Ai_{kn} <= 1, i=1,...,m,
%              hi*x1^G_{i1}*...*xn^G_{in} = 1,            i=1,...,p,
%              li <= xi <= ui,                            i=1,...,n,
%
% where variables are x1,...,xn and the problem data are bi_k, Ai_{kj},
% for i = 1,...,m, hi, G_{ij} for i = 1,...,p.
%
% Calling sequences:
%
%  [x,status,lambda,nu] = gpposy(A,b,szs)
%  [x,status,lambda,nu] = gpposy(A,b,szs,G,h)
%  [x,status,lambda,nu] = gpposy(A,b,szs,G,h,l,u)
%  [x,status,lambda,nu] = gpposy(A,b,szs,G,h,l,u,quiet)
%
% Examples:
%
%  [x,status,lambda,nu] = gpposy(A,b,szs,G,h)
%  [x,status,lambda,nu] = gpposy(A,b,szs,[],[],[],[],quiet)
%
% Input arguments:
%
% - A:         (sum_i n_i) x n matrix;  A   = [A0' A1' ... Am' ]'
% - b:         (sum_i n_i) vector;      b   = [b0' b1' ... bm' ]'
% - szs:       dimensions of Ai and bi; szs = [n0 n1 ... nm]'
%              where Ai is (ni x n) and bi is (ni x 1)
% - G:         p x n matrix
% - h:         p-vector
% - l:         n-vector; lower bound for x
% - u:         n-vector; upper bound for x
% - quiet:     boolean variable; suppress all the print messages if true
%
% Output arguments:
%
% - x:         n-vector; primal optimal point
% - nu:        (sum_i n_i) vector;  nu = [nu0' nu1' ... num']'
%              optimal sensitivity for constraints Ai*x + bi = yi
% - mu:        (sum_i l_i) vector;  mu = [mu0' mu1' ... mul']'
%              optimal sensitivity for constraints G*x + h = 0
% - lambda:    m-vector, optimal sensitivity for inequality constraints
% - status:    string; 'INFEASIBLE', 'SOLVED', or 'FAILED'
%              
% gpposy changes the problem data into posynomial form and calls the the
% sovler, gpcvx.m.

%----------------------------------------------------------------------
%       INITIALIZATION
%----------------------------------------------------------------------

% PARAMETERS
LOWER_BOUND = 1e-100;
UPPER_BOUND = 1e+100;

n   = size(A,2);     % # of variables (x1,..,xn)

% VARIABLE ARGUMENT HANDLING
defaults  = {[],[],LOWER_BOUND*ones(n,1),UPPER_BOUND*ones(n,1),false};
givenArgs = ~cellfun('isempty',varargin);
defaults(givenArgs) = varargin(givenArgs);
[G,h,l,u,quiet]  = deal(defaults{:});

% CONVERT PROBLEM DATA INTO CONVEX FORM
b = log(b); h = log(h); l = log(l); u = log(u);

if ( ~(all(isfinite(b)) && all(isfinite(h))) )
    disp('ERROR: Too small value of b or h');
    x = []; status = 'FAILED'; lambda = []; nu = [];
    return;
end
%----------------------------------------------------------------------
%       Call gpcvx
%----------------------------------------------------------------------
[x,status,lambda,empty,nu] = gpcvx(A,b,szs,G,h,l,u,quiet);

% CONVERT SOLUTION TO POSYNOMIAL FORM
x = exp(x);

%**********************************************************************
% helper functions for gpposy.m
% (merged by almir into a single file, 12/6/06)
%**********************************************************************

%**********************************************************************
% gpcvx.m
%**********************************************************************
function [x,status,lambda,nu,mu] = gpcvx(A,b,szs,varargin)

% [x,status,lambda,nu,mu] = gpcvx(A,b,szs,G,h,l,u,quiet)
%
% solves the geometric program in convex form
%
%  minimize    lse(y0)
%  subject to  lse(yi) <= 0,   i=1,...,m,
%              Ai*x+bi = yi,   i=0,...,m,
%              G*x+h = 0,
%              li <= xi <= ui, i=1,...,n
%
% where lse is defined as  lse(y) = log sum_i exp yi,
% and the dual problem,
%
%  maximize    b0'*nu0 + ... + bm'*num + h'*mu + lambdal'*l - lambdau'*u +
%                 entr(nu0) + lambda1*entr(nu1/lambda1) + 
%                 ,..., + lambdam*entr(num/lambdam)
%  subject to  nui >= 0,         i=0,...,m,
%              lambdai >= 0,     i=1,...,m,
%              1'*nu0 = 1
%              1'*nui = lambdai, i=1,...,m,
%              A0'*nu0 + ... + Am'*num + G'*mu + lambdau - lambdal = 0,
%
% where entr is defined as  entr(y) = -sum_i yi*log(yi).
%
% Calling sequences:
%
%  [x,status,lambda,nu,mu] = gpcvx(A,b,szs)
%  [x,status,lambda,nu,mu] = gpcvx(A,b,szs,G,h)
%  [x,status,lambda,nu,mu] = gpcvx(A,b,szs,G,h,l,u)
%  [x,status,lambda,nu,mu] = gpcvx(A,b,szs,G,h,l,u,quiet)
%
% Examples:
%
%  [x,status,lambda,nu,mu] = gpcvx(A,b,szs,G,h)
%  [x,status,lambda,nu,mu] = gpcvx(A,b,szs,[],[],[],[],quite)
%
% Input arguments:
%
% - A:         (sum_i n_i) x n matrix;  A   = [A0' A1' ... Am' ]'
% - b:         (sum_i n_i) vector;      b   = [b0' b1' ... bm' ]'
% - szs:       dimensions of Ai and bi; szs = [n0 n1 ... nm]'
%              where Ai is (ni x n) and bi is (ni x 1)
% - G:         p x n matrix
% - h:         p-vector
% - l:         n-vector; lower bound for x
% - u:         n-vector; upper bound for x
% - quiet:     boolean variable; suppress all the print messages if true
%
% Output arguments:
%
% - x:         n-vector; primal optimal point
% - nu:        (sum_i n_i) vector;  nu = [nu0' nu1' ... num']'
%              dual variables for constraints Ai*x + bi = yi
% - mu:        (sum_i l_i) vector;  mu = [mu0' mu1' ... mul']'
%              dual variables for constraints G*x + h = 0
% - lambda:    m-vector, dual variables for inequality constraints
% - status:    string; 'Infeasible', 'Solved', or 'Failed'
%              
% gpcvx sets up phase 1 and phase 2 and calls the real sovler, gppd2.m.

%----------------------------------------------------------------------
%       INITIALIZATION
%----------------------------------------------------------------------

% PARAMETERS
LOWER_BOUND = -250;
UPPER_BOUND = +250;

% DIMENSIONS
N  = size(A,1);     % # of terms in the obj. and inequalities
n  = size(A,2);     % # of variables (x1,..,xn)
m  = length(szs)-1; % # of inequalities
n0 = szs(1);        % # of terms in the objective

% VARIABLE ARGUMENT HANDLING
defaults  = {[],[],LOWER_BOUND*ones(n,1),UPPER_BOUND*ones(n,1),false};
givenArgs = ~cellfun('isempty',varargin);
defaults(givenArgs) = varargin(givenArgs);
[G,h,l,u,quiet]  = deal(defaults{:});

% MATLAB LIMIT OF LOWER/UPPER BOUNDS
l(l<LOWER_BOUND) = LOWER_BOUND;
u(u>UPPER_BOUND) = UPPER_BOUND;

if (isempty(G)) 
    G = zeros(0,n);
    h = zeros(0,1);
end
p  = size(G,1);     % # of (terms in) the equality constraints

% E is a matrix s.t. [1'*y0  1'*y1  ... 1'*ym ]' = E*y
indsl = cumsum(szs);
indsf = indsl-szs+1;
lx    = zeros(N,1);
lx(indsf) = 1;
E = sparse(cumsum(lx),[1:N],ones(N,1));

%----------------------------------------------------------------------
%               PHASE I
%----------------------------------------------------------------------

% solves the feasibility problem
%
%  minimize    s
%  subject to  lse(yi) <= s,     i=1,...,m,
%              Ai*x+bi = yi,     i=0,...,m,
%              G*x+h = 0,
%              li <= xi <= ui,   i=1,...,n
%
% where lse is defined as  lse(y) = log sum_i exp yi,
%
% For phase I
% 1) change objective function to s
% 2) change constraints from fi(x) <= 0 to fi(x) <= s
% 3) add bound constraints; li <= xi <= ui
%
% Hence, we set up a new objective and constraints, 
%    i.e., A,b,G,h and szs for Phase I optimization.
%
% Change the size vector, szs.
%
% Change A and b
%          s    xi
%   A1 = [ 1 | 0 0 0      b1 = [ 0      <- (a1) new objective
%         ---+------             -
%         -1 |
%         -1 |A(ineq)            b      <- (a2) new inequalities
%         -1 |
%         ---+------             -
%         -1 |-1 0 0             l1
%         -1 | 0-1 0             l2     <- (a4) l <= xi
%         -1 | 0 0-1             l3
%         ---+------             -
%         -1 | 1 0 0            -u1
%         -1 | 0 1 0            -u2     <- (a6) xi <= u
%         -1 | 0 0 1 ];         -u3 ];

% FORM SZS
szs1    = [1; szs(2:end); ones(2*n,1) ];

% FORM INITIAL X
if (p == 0)
    xinit = zeros(n,1);
else
    xinit = G'*((G*G')\h);
end

% FORM INITIAL S
%  sinit = max(fi,0) since fi <= si and 0 <= si

y = A*xinit+b;
[f,expyy] = lse(E,y);
finit = f(2:m+1);
linit = -xinit+l;
uinit = +xinit-u;
sinit = max([0; finit; linit; uinit]) + 1; % + 1 is for margin.

% FORM A AND B
A1 = [+1            , sparse(1,n);...
      -ones(N-n0,1) , A(n0+1:N,:);...
      -ones(n,1)    ,-speye(n)   ;...
      -ones(n,1)    ,+speye(n)   ];
b1 = [ 0; b((n0+1):N); l; -u ];

% FORM G AND H
G1 = [spalloc(size(G,1),1,0), G];
h1 = [h];

% CALL THE INTERNAL GP SOLVER
[x,status,lambda,nu,mu] = gppd2(A1,b1,szs1,[sinit;xinit],G1,h1,true,quiet);

% EXTRACT X FROM [S; X]
x0 = x(2:n+1);

y = A*x0+b;
[f,expyy] = lse(E,y);
f1m       = f(2:m+1);

% FEASIBILITY CHECK OF PHASE I SOLUTION
if (status <= 0 || max([f1m; -Inf]) >= 0)
    status = 'Infeasible';
    if (~quiet) disp(status); end
    return
end
clear A1 b1 G1 h1 x01 szs1;        

%----------------------------------------------------------------------
%               PHASE II
%----------------------------------------------------------------------

% solves the geometric program in convex form
%
%  minimize    lse(y0)
%  subject to  lse(yi) <= 0,   i=1,...,m,
%              Ai*x+bi = yi,   i=0,...,m,
%              G*x+h = 0,
%              li <= xi <= ui, i=1,...,n
%
% where lse is defined as  lse(y) = log sum_i exp yi,
%
% Change A and b to add the bound of x into the inequality constraints.
%           
%   A2 = [  A         b2 = [ b
%         ------            ----
%         -1 0 0             l1
%          0-1 0             l2     <- li <= xi
%          0 0-1             l3
%         ------            ----
%          1 0 0            -u1     <- xi <= ui
%          0 1 0            -u2
%          0 0 1 ];         -u3 ];

szs2 = [ szs ; ones(2*n,1) ];
A2   = [ sparse(A); -speye(n); speye(n) ];
b2   = [ b; l; -u ];

% CALL THE INTERNAL GP SOLVER
[x,status,lambda,nu,mu] = gppd2(A2,b2,szs2,x0,G,h,false,quiet);

if (status <= 0)
    status = 'Failed';
    if (~quiet) disp(status); end
    return
else
    status = 'Solved';
    if (~quiet) disp(status); end
    return
end

%**********************************************************************
% gppd2.m
%**********************************************************************
function [x,status,la,nu,mu] = gppd2(A,b,szs,x0,G,h,phase1,quiet)

% [x,nu,mu,la,status] = gppd2(A,b,szs,x0,G,h,phase1,quiet)
%
% solves the geometric program in convex form with a starting point.
%
%  minimize    lse(y0)
%  subject to  lse(yi) <= 0,   i=1,...,m,
%              Ai*x+bi = yi,   i=0,...,m,
%              G*x+h = 0,
%
% where lse is defined as  lse(y) = log sum_i exp yi,
%
% and the dual problem,
%
%  maximize    b0'*nu0 + ... + bm'*num + h'*mu +
%                 entr(nu0) + la1*entr(nu1/la1) + 
%                 ,..., + lam*entr(num/lam)
%  subject to  nui >= 0,         i=0,...,m,
%              lai >= 0,     i=1,...,m,
%              1'*nu0 = 1
%              1'*nui = lai, i=1,...,m,
%              A0'*nu0 + ... + Am'*num + G'*mu = 0,
%
% where entr is defined as  entr(y) = -sum_i yi*log(yi).
%
% x0 should satisfy the primal inequality constraints.
%
% Input arguments:
%
% - A:         (sum_i n_i) x n matrix; A = [A0' A1' ... Am' ]'
% - b:         (sum_i n_i) vector;   b = [b0' b1' ... bm' ]'
% - szs:       dimensions of Ai and bi; szs = [Nob n1 ... nm]' 
%              where Ai is (ni x n) and bi is (ni x 1)
% - x0:        n-vector; MUST BE STRICTLY FEASIBLE FOR INEQUALITIES
% - G:         p x n matrix
% - h:         p-vector
% - phase1:    boolean variable; indicator for phase I and phase II
%              true -> Phase I, false -> Phase II
% - quiet:     boolean variable; suppress all the print messages if true
%
% Output arguments:
%
% - x:         n-vector; primal optimal point
% - nu:        (sum_i n_i) vector;  nu = [nu0' nu1' ... num']'
%              dual variables for constraints Ai*x + bi = yi
% - mu:        p vector; mu = [mu1 ... mup]'
%              dual variables for constraints G*x + h = 0
% - la:        m vector; la = [lambda1 ... lambdam]'
%              dual variables lambda; la_i = sum(nu_i)
% - status:    scalar;
%              2	Function converged to a solution x.
%              1	Phase I success; x(1:szs(1)) <= 0.
%             -1	Number of iterations exceeded MAXITERS.
%             -2	Starting point is not strictly feasible.
%             -3	Newton step calculation failure.

%----------------------------------------------------------------------
%               INITIALIZATION
%----------------------------------------------------------------------

% PARAMETERS
ALPHA   = 0.01;     % backtracking linesearch parameter (0,0.5]
BETA    = 0.5;      % backtracking linesearch parameter (0,1)
MU      = 2;        % IPM parameter: t update
MAXITER = 500;      % IPM parameter: max iteration of IPM
EPS     = 1e-5;     % IPM parameter: tolerance of surrogate duality gap
EPSfeas = 1e-5;     % IPM parameter: tolerance of feasibility

% DIMENSIONS
[N,n] = size(A); m = length(szs)-1; p = size(G,1); n0 = szs(1);
if (isempty(G)), G = zeros(0,n); h = zeros(0,1); end

warning off all;

% SPARSE ZERO MATRIX
Opxp = sparse(p,p);

% SUM MATRIX: E is a matrix s.t. [1'*y0  1'*y1  ... 1'*ym ]' = E*y
indsl = cumsum(szs); 
indsf = indsl-szs+1;
lx    = zeros(N,1);
lx(indsf) = 1;
E = sparse(cumsum(lx),[1:N],ones(N,1));

x = x0;
% f1m is a LHS vector of inequality constraints s.t [f1' ... fm']'
y = A*x+b;
[f,expyy] = lse(E,y);
f1m       = f(2:m+1);

% CHECK THE INITIAL CONDITIONS
if (max(f1m) >= 0)
   if (~quiet) disp(['x0 is not strictly feasible.']); end
   la = []; nu = []; mu = [];
   status = -2;
   return;
end;

% INITIAL DUAL POINT
la = -1./f1m;   % positive value with duality gap 1.
nu = ones(N,1); % ANY value.
mu = ones(p,1); % ANY value.

step = Inf;
pp = [];

if (~quiet) disp(sprintf('\n%s %15s %11s %20s %18s \n',...
    'Iteration','primal obj.','gap','dual residual','previous step.')); end

%----------------------------------------------------------------------
%               MAIN LOOP
%----------------------------------------------------------------------
for iters = 1:MAXITER

    gap = -f1m'*la;

    % UPDATE T
    % update t only when the current x is not to far from the cetural path.
    if (step > 0.5)
        t = m*MU/gap;
    end

    % CALCULATE RESIDUALS
    % gradfy = exp(y)./(E'*(E*exp(y)));
    gradfy   = expyy./(E'*(E*expyy));
    resDual  = A'*(gradfy.*(E'*[1;la])) + G'*mu;
    resPrim  = G*x + h;
    resCent  = [-la.*f1m-1/t];
    residual = [resDual; resCent; resPrim];

    if (~quiet) disp(sprintf('%4d %20.5e %16.5e %14.2e %16.2e',...
        iters,f(1),-f1m'*la,norm(resDual),step)); end

    % STOPPING CRITERION FOR PHASE I
    if ((phase1 == true) & max(x(1:szs(1))) < 0)
        nu = gradfy.*(E'*[1;la]);
        status = 1;
        return;
    end;
    % STOPPING CRITERION FOR PHASE I & II
    if ( (gap <= EPS) & ...
            (norm(resDual) <= EPSfeas) & (norm(resPrim) <= EPSfeas) )
        % this calculation of nu is correct only when reached to optimal
        nu = gradfy.*(E'*[1;la]);
        status = 2;
        return;
    end;

    % CALCULATE NEWTON STEP
    diagL1  = sparse(1:m+1,1:m+1,[0;-la./f1m]);
    diagL2  = sparse(1:m+1,1:m+1,[1;la]);
    diagG1  = sparse(1:N,1:N,gradfy);
    diagG2  = sparse(1:N,1:N,gradfy.*(E'*[1;la]));

    EGA = E*diagG1*A;
    H1  = EGA'*(diagL1-diagL2)*EGA + A'*diagG2*A;

	if( 0 )
        dz = -[ H1, G'   ; ...
                G , Opxp ] ...
               \ ...
               [EGA'*[1;-1./(t*f1m)]+G'*mu ; resPrim ];
    else
        dz = -H1 \ [EGA'*[1;-1./(t*f1m)]];
	end

    % PERTURB KKT WHEN (ALMOST) SINGULAR
    %  add small diagonal terms when the matrix is almost singular.
    perturb = EPS;
    while (any(isinf(dz)) || any(isnan(dz)))
       dz = -([ H1, G'; G , Opxp ] + sparse(1:n+p,1:n+p,perturb))...
             \ ...
             [EGA'*[1;-1./(t*f1m)]+G'*mu ; resPrim ];
        % increase the size of diagonal term if still singular
        perturb = perturb*10;
        if (perturb > 1) break; end
    end
    % ERROR CHECK FOR NEWTON STEP CALCULATION
    if (any(isnan(dz)))
        nu = gradfy.*(E'*[1;la]);
        status = -3;
        return;
    end

    dx  = dz(1:n); dy  = A*dx; dmu = dz(n+[1:p]');
    dla = -la./f1m.*(E(2:m+1,:)*(gradfy.*dy))+resCent./f1m;

    % BACKTRACKING LINESEARCH
    negIdx = (dla < 0);
    if (any(negIdx))
        step = min( 1, 0.99*min(-la(negIdx)./dla(negIdx)) );
    else
        step = 1;
    end
    while (1)
        newx    = x  + step*dx;
        newy    = y  + step*dy;
        newla   = la + step*dla;
        newmu   = mu + step*dmu;
        [newf,newexpyy] = lse(E,newy);
        newf1m  = newf(2:end);

        % UPDATE RESIDUAL
        % newGradfy = exp(newy)./(E'*(E*exp(newy)));
        newGradfy = newexpyy./(E'*(E*newexpyy));

        newResDual  = A'*(newGradfy.*(E'*[1;newla])) + G'*newmu;
        newResPrim  = G*newx + h;
        newResCent  = [-newla.*newf1m-1/t];
        newResidual = [newResDual; newResCent; newResPrim];
        
        if ( max(newf1m) < 0 && ...
             norm(newResidual) <= (1-ALPHA*step)*norm(residual) )
            break;
        end
        step = BETA*step;
    end;
    % UPDATE PRIMAL AND DUAL VARIABLES
    x  = newx; mu = newmu; la = newla; y = A*x+b;
    [f,expyy] = lse(E,y); f1m = f(2:m+1);
end
if (iters >= MAXITER)
    if (~quiet) disp(['Maxiters exceeded.']); end
    nu = gradfy.*(E'*[1;la]);
    status = -1;
    return;
end

%**********************************************************************
% lse.m
%**********************************************************************
function [f,expyy] = lse(E,y)
% 
% LSE finds log-sum-exp values.
%
% It calculates log(E*exp(y)), with special handling of extreme y values.
% 

% ymax is a vector of maximum exponent values of each posynomial.
%
% ex: log(exp(y1)+exp(y2)+exp(y3)), log(exp(y4)+exp(y5))
%     ymax = [max(y1,y2,y3); max(y4,y5)]
%
% note: min(y) are substracted and added to handle the case when all the
%       exponents of a posynomial are negative.
%       i.e., 0 0 0 -3 -4 -> max value should be -3

ymax  = full(max(E*sparse(1:length(y),1:length(y),y-min(y)),[],2))+min(y);
expyy = exp(y-E'*ymax);
f     = log(E*expyy)+ymax;
