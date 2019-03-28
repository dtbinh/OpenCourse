function [s,lams,iters,gaps,residl] = pd_thr_speed(G,phi,Tother,Tamb,Tmax,smin,smax,TOL,quiet,varargin)
% [s,lam,iters,gaps,resdls] = pd_thr_speed(G,phi,Tother,Tamb,Tmax,smin,smax,TOL,quiet,varargin)
%
% This is implementation of the infeasible-start primal-dual interior-point
% method described in "Processor speed control with thermal constraints"
% by Mutapcic, Boyd, Murali, Atienza, De Micheli, Gupta
%
% Algorithm solves the optimal processor speed scaling with thermal constraints:
%
%    maximize    U(s)
%    subject to  G\phi(s) + Tother + Tamb <= Tmax
%                s_min <= s <= s_max
%
% where processor speeds s in R^n are the optimization variables,
% and the problem data is listed below. For more detail, consult the paper.
%
% Input arguments:
%
% - G:      thermal marix (can pass a sparse matrix in too)
% - phi:    functions handles that give \phi, grad(\phi), and Hess(\phi), e.g.,
%           phi.fval = @(s) s.^3;
%           phi.grad = @(s) 3*s.^2;
%           phi.hess = @(s) spdiags(6*s,0,length(s),length(s));
% - Tother: other temperatures (vector in R^m)
% - Tamb:   ambient temperature (a real)
% - Tmax:   maximum allowable temperature (a real)
% - smin:   minimum speed (a real)
% - smax:   maximum speed (a real)
% - TOL:    tolerance for the duality measure
% - quiet:  be quiet (set to 1) or verbose and report results (set to 0)
%
% where varargin can include:
%
% - s0:     initial feasible processor speed (vector in R^n)
% - lams0:  initial feasible dual variables, all of them (vector in R^(m+2n), see paper)

% Almir Mutapcic, 9/07

[m,n] = size(G); N = m+2*n;

% SLATER'S CONDITION CHECK
z  = Tmax - G*phi.fval(smin*ones(n,1)) - Tother - Tamb;
if any(z <= 0), error('Slater condition is violated. Exiting.'); end

% CONSTANTS
MAXITERS = 40;
GAMMA    = 0.05;
BETA     = 5;
SIGMAMIN = 0.15;
SIGMAMAX = 0.5;
REDUCE   = 0.85;

Inxn = speye(n,n);
Imxm = speye(m,m);

Onxn = sparse(n,n);
Omxm = sparse(m,m);
Omxn = sparse(m,n);
Onxm = sparse(n,m);

% INITIAL VALUES
if ( nargin == 11 )
  [s lams] = deal( varargin{:} );
  lam  = lams(1:m);
  lamu = lams(m+1:m+n);
  laml = lams(m+n+1:end);
elseif ( nargin == 9 )
  s   = (smin + 1e-3)*ones(n,1);
else
  error('Wrong number of arguments: type help pd_thr_speed')
end

z  = Tmax - G*phi.fval(s) - Tother - Tamb;
zu = smax - s;
zl = s - smin;

while (min([z;zu;zl]) <= 0), 
   disp(['s0 is not strictly feasible.']);
   s  = max( smin+1e-5, s*0.995 ); 
   z  = Tmax - G*phi.fval(s) - Tother -Tamb;
   zu = smax - s;
   zl = s - smin;
end;

if ( nargin == 9 )
  % start with gap m+2*n
  lam  = 1./z;
  lamu = 1./zu;
  laml = 1./zl;
end


gaps = []; residl = [];
optval = sum(s);
mu = (lam'*z + lamu'*zu + laml'*zl)/N;  gap = mu;

resPri  = G*phi.fval(s) + Tother + Tamb - Tmax + z;
resDual = ones(n,1) - spdiags(phi.grad(s),0,n,n)*G'*lam - lamu + laml;
ratioDual = norm(resDual)/mu;


step = Inf;
if ~quiet
  disp([' ']);
  fprintf(1,'Iteration     primal obj.          gap         primal res.');
  fprintf(1,'      dual res.        step.\n');
end


for iters = 1:MAXITERS


   if ~quiet
      disp([sprintf('%4d',iters),  sprintf('% 20.5e',optval),  ...
            sprintf('% 18.5e', gap),  sprintf('% 15.2e', norm(resPri)), ...
            sprintf('% 15.2e', norm(resDual)), sprintf('% 18.5e', step) ]);
   end


   % SET MU AND SIGMA
   mu = (lam'*z + lamu'*zu + laml'*zl)/N;  gap = mu;
   sigma = SIGMAMAX;


   % STOPPING CRITERION
   gaps = [gaps gap];
   residl = [residl norm([resPri; resDual])];
   if ( gap < TOL )
      lams = [lam; lamu; laml];
      return;
   end;


   % SOLVING NEWTON SYSTEM USING BACK-SUBSTITUTION
   D     = spdiags(phi.grad(s),0,n,n);
   DGT   = D*G';
   Hs    = phi.hess(s)*spdiags(G'*lam,0,n,n);

   resDual  = ones(n,1) - DGT*lam - lamu + laml;
   resPri   = G*phi.fval(s) + Tother + Tamb - Tmax + z;
   resU     = s - smax + zu;
   resL     = smin - s + zl;
   resCentP = [ lam.*z   ] - sigma*mu;
   resCentU = [ lamu.*zu ] - sigma*mu;
   resCentL = [ laml.*zl ] - sigma*mu;

   % first find optimal ds
   H  = Hs + DGT*spdiags(lam./z,0,m,m)*DGT' + spdiags(lamu./zu,0,n,n) + spdiags(laml./zl,0,n,n);
   r  = resDual + DGT*(resCentP./z - lam./z.*resPri) + ...
        resCentU./zu - lamu./zu.*resU - resCentL./zl + laml./zl.*resL;

   % COMPUTE THE NEWTON STEP (BACKSLASH USES CHOLESKY HERE)
   ds = H \ r;

   dz  = -DGT'*ds - resPri;
   dzu = -ds - resU;
   dzl = +ds - resL;

   % solve for lambda dual variables
   dlam  = lam./z.*(DGT'*ds + resPri) - resCentP./z;
   dlamu = lamu./zu.*(+ds + resU) - resCentU./zu;
   dlaml = laml./zl.*(-ds + resL) - resCentL./zl;


   % LINE SEARCH
   neglam  = (dlam  < 0);
   neglamu = (dlamu < 0);
   neglaml = (dlaml < 0);
   negz    = (dz  < 0);
   negzu   = (dzu < 0);
   negzl   = (dzl < 0);

   step = 1;
   if (any(neglam))
      step = min( step, 0.99*min(-lam(neglam)./dlam(neglam)) );
   end
   if (any(neglamu))
      step = min( step, 0.99*min(-lamu(neglamu)./dlamu(neglamu)) );
   end
   if (any(neglaml))
      step = min( step, 0.99*min(-laml(neglaml)./dlaml(neglaml)) );
   end
   if (any(negz))
      step = min( step, 0.99*min(-z(negz)./dz(negz)) );
   end
   if (any(negzu))
      step = min( step, 0.99*min(-zu(negzu)./dzu(negzu)) );
   end
   if (any(negzl))
      step = min( step, 0.99*min(-zl(negzl)./dzl(negzl)) );
   end


   news    = s    + step*ds;
   newlam  = lam  + step*dlam;
   newlamu = lamu + step*dlamu;
   newlaml = laml + step*dlaml;
   newz    = z    + step*dz;
   newzu   = zu   + step*dzu;
   newzl   = zl   + step*dzl;

   newmu      = (newlam'*newz + newlamu'*newzu + newlaml'*newzl)/N;
   newresDual = ones(n,1) - spdiags(phi.grad(news),0,n,n)*G'*newlam - newlamu + newlaml;
   newresPri  = G*phi.fval(news) + Tother + Tamb - Tmax + newz;
   newresU    = news - smax + newzu;
   newresL    = smin - news + newzl;

   while( min( [newlam.*newz; newlamu.*newzu; newlaml.*newzl] ) < GAMMA*newmu & ...
          norm(newresDual) > ratioDual*BETA*newmu & ...
		  norm([newresPri; newresU; newresL;]) > BETA*newmu & ...
          newmu > (1-0.01*step)*mu )

      step    = REDUCE*step;
      news    = s    + step*ds;
      newlam  = lam  + step*dlam;
      newlamu = lamu + step*dlamu;
      newlaml = laml + step*dlaml;
      newz    = z    + step*dz;
      newzu   = zu   + step*dzu;
      newzl   = zl   + step*dzl;

      newmu      = (newlam'*newz + newlamu'*newzu + newlaml'*newzl)/N;
      newresDual = ones(n,1) - spdiags(phi.grad(news),0,n,n)*G'*newlam - newlamu + newlaml;
      newresPri  = G*phi.fval(news) + Tother + Tamb - Tmax + newz;
      newresU    = news - smax + newzu;
      newresL    = smin - news + newzl;
   end

   % UPDATE THE SOLUTION
   s    = s    + step*ds;
   lam  = lam  + step*dlam;
   lamu = lamu + step*dlamu;
   laml = laml + step*dlaml;
   z    = z    + step*dz;
   zu   = zu   + step*dzu;
   zl   = zl   + step*dzl;

   % exact computation of slacks z (also works well)
   %% z    = Tmax - G*phi.fval(s) - Tother - Tamb;
   %% zu   = smax - s;
   %% zl   = s - smin;

   mu      = (lam'*z + lamu'*zu + laml'*zl)/N;  gap = mu;
   resDual = ones(n,1) - spdiags(phi.grad(s),0,n,n)*G'*lam - lamu + laml;
   optval = sum(s);
 
end;

% FINAL CHECK IF MAXITERS EXCEEDED
if (iters == MAXITERS),
   disp(['Maxiters exceeded.']);
   s=[];  lams=[];  iters=[];  gaps=[];  residl=[];
   return;
end;
