function [x,status] = l1tf(y,lambda)
% [x,status] = l1tf(y,lambda)
%
% finds the solution of the l1 trend estimation problem
%
%  minimize    (1/2)||y-x||^2+lambda*||Dx||_1,
%
% with variable x, and problem data y and lambda, with lambda >0.
% D is the second difference matrix, with rows [0... -1 2 -1 ...0]
%
% and the dual problem
%
%  minimize    (1/2)||D'*z||^2-y'*D'*z
%  subject to  norm(z,inf) <= lambda,
%
% with variable z.
%
% Input arguments:
%
% - y:          n-vector; original signal
% - lambda:     scalar; positive regularization parameter
%
% Output arguments:
%
% - x:          n-vector; primal optimal point
% - z:          2(n-1)-vector; dual optimal point
% - status:     string;
%               'solved', 'maxiter exceeded'
%
% for more details,
% see "l1 Trend Filtering", S. Kim, K. Koh, ,S. Boyd and D. Gorinevsky
% www.stanford.edu/~boyd/l1_trend_filtering.html
%
 
%----------------------------------------------------------------------
%               INITIALIZATION
%----------------------------------------------------------------------

% PARAMETERS
ALPHA     = 0.01;   % backtracking linesearch parameter (0,0.5]
BETA      = 0.5;    % backtracking linesearch parameter (0,1)
MU        = 2;      % IPM parameter: t update
MAXITER   = 40;     % IPM parameter: max iteration of IPM
MAXLSITER = 20;     % IPM parameter: max iteration of line search
TOL       = 1e-4;   % IPM parameter: tolerance

% DIMENSIONS
n   = length(y);    % length of signal x
m   = n-2;          % length of Dx

%warning off all;

% OPERATOR MATRICES
I2  = speye(n-2,n-2);
O2  = zeros(n-2,1);
D   = [I2 O2 O2]+[O2 -2*I2 O2]+[O2 O2 I2];

DDT = D*D';
Dy  = D*y;

% VARIABLES
z   = zeros(m,1);   % dual variable
mu1 = ones(m,1);    % dual of dual variable
mu2 = ones(m,1);    % dual of dual variable

t    = 1e-10; 
pobj =  Inf;
dobj =  0;
step =  Inf;
f1   =  z-lambda;
f2   = -z-lambda;

disp('--------------------------------------------');
disp('l1 trend filtering via primal-dual algorithm');
disp('Version 0.7 May 1 2007');
disp('Kwangmoo Koh, Seung-Jean Kim, Stephen Boyd');
disp('--------------------------------------------');
%disp(sprintf('\n%s %13s %12s %8s %9s %17s \n',...
    %'Iteration','Primal obj.','Dual obj.','Gap','t','Step size'));
disp(sprintf('\n%s %13s %12s %8s\n',...
    'Iteration','Primal obj.','Dual obj.','Gap'));

%----------------------------------------------------------------------
%               MAIN LOOP
%----------------------------------------------------------------------
for iters = 0:MAXITER

    DTz  = (z'*D)';
    DDTz = D*DTz;
    w    = Dy-(mu1-mu2);
    
    % two ways to evaluate primal objective:
    % 1) using dual variable of dual problem
    % 2) using optimality condition
    pobj1 = 0.5*w'*(DDT\w)+lambda*sum(mu1+mu2);
    pobj2 = 0.5*DTz'*DTz+lambda*sum(abs(Dy-DDTz));
    pobj = min(pobj1,pobj2);
    dobj = -0.5*DTz'*DTz+Dy'*z;
%     pobj = min(min(pobj1,pobj2), pobj);
%     dobj = max(-0.5*DTz'*DTz+Dy'*z, dobj);
    gap  =  pobj - dobj;

    %disp(sprintf('%6d %16.4e %13.5e %10.2e %11.2e %13.2e',...
        %iters, pobj, dobj, gap, t, step));
    disp(sprintf('%6d %15.4e %13.5e %10.2e',...
        iters, pobj, dobj, gap));

    % STOPPING CRITERION
    if (gap <= TOL)
        status = 'solved';
        disp(status);
        x = y-D'*z;
        return;
    end;

    if (step >= 0.2)
        t =max(2*m*MU/gap, 1.2*t);
    end

    % CALCULATE NEWTON STEP
    
    rz      =  DDTz - w;
    S       =  DDT-sparse(1:m,1:m,mu1./f1+mu2./f2);
    r       = -DDTz + Dy + (1/t)./f1 - (1/t)./f2;
    dz      =  S\r;
    dmu1    = -(mu1+((1/t)+dz.*mu1)./f1);
    dmu2    = -(mu2+((1/t)-dz.*mu2)./f2);

    resDual = rz;
    resCent = [-mu1.*f1-1/t; -mu2.*f2-1/t];
    residual= [resDual; resCent];

    % BACKTRACKING LINESEARCH
    negIdx1 = (dmu1 < 0); 
    negIdx2 = (dmu2 < 0);
    step = 1;
    if (any(negIdx1))
        step = min( step, 0.99*min(-mu1(negIdx1)./dmu1(negIdx1)) );
    end
    if (any(negIdx2))
        step = min( step, 0.99*min(-mu2(negIdx2)./dmu2(negIdx2)) );
    end

    for liter = 1:MAXLSITER
        newz    =  z  + step*dz;
        newmu1  =  mu1 + step*dmu1;
        newmu2  =  mu2 + step*dmu2;
        newf1   =  newz - lambda;
        newf2   = -newz - lambda;

        % UPDATE RESIDUAL
        newResDual  = DDT*newz - Dy + newmu1 - newmu2;
        newResCent  = [-newmu1.*newf1-1/t; -newmu2.*newf2-1/t];
        newResidual = [newResDual; newResCent];
        
        if ( max(max(newf1),max(newf2)) < 0 && ...
            norm(newResidual) <= (1-ALPHA*step)*norm(residual) )
            break;
        end
        step = BETA*step;
    end
    % UPDATE PRIMAL AND DUAL VARIABLES
    z  = newz; mu1 = newmu1; mu2 = newmu2; f1 = newf1; f2 = newf2;
end

% The solution may be close at this point, but does not meet the stopping
% criterion (in terms of duality gap).
x = y-D'*z;
if (iters >= MAXITER)
    status = 'maxiter exceeded';
    disp(status);
    return;
end
