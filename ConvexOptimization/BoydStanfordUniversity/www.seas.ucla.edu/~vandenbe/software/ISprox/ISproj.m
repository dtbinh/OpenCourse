function [sol,gphi,phi,stats] = ISproj(u,leq,DC,tol)
%  (The procedure calls the function toeplitzFact)
%  [sol,gphi,phi,stats] = ISPROJ(u,leq,DC,tol)
%  Solves the following problem when leq == false
%  min.  -<u,x> + phi(x)
%  s.t.  x_0 = DC.
%  Solves the problem with equality constraint replaced by
%  x_0 <= DC when leq == true.
%  The default DC is 1 when not specified.
%  The default tol is 1e-8 when not specified.
%  
%  Note:  <u,x> = u0 x0 + 2 u1 x1 + 2 u2 x2 + ... + 2 up xp
%         X(w) = x_0 + 2 sum_{k=1}^p x_k cos(kw)
%         phi(x) = - (1/2 pi) \int_{-pi}^{pi} log X(w) dw
%
%  INPUT:
%         u:  real vector
%       leq:  boolean 
%        DC:  positive scalar
%       tol:  tolerance for Newton's method, default tol = 1e-8
%  OUTPUT:
%       sol:  the solution x
%      gphi:  gradient of phi at x
%       phi:  function value phi(x) at x
%     stats:  procedure statistics
%
%  Author: Hsiao-Han (Suzi) Chao (suzihchao@ucla.edu)
%  Date modified: 7/18/2018

stats = zeros(4,1); %[#iter; #newton; #levinson; cputime]
stats(4) = cputime;

DEBUG = 0;

if nargin < 4
    tol = 1e-8; 
end
if nargin == 2
    DC = 1;
end
p = length(u) - 1;
maxiter = 50;

%% Safeguarded Newton's method
% use toeplitzFact to check positive definiteness and get a decomposition
% and use a hybrid of bisection and Newtons method to find v | |f'(v)| < tol
% for vI - T(u) > 0
% f(v) = -log(e0'(vI-T(u)^(-1)e0) - DC*v
% f'(v) = e0'(vI-T(u))^(-2)e0/e0'(vI-T(u))^(-1)e0 - DC
% f''(v) = ((e0'(vI-T(u))^(-2)e0)^2 -
%          2*e0'(vI-T(u))^(-1)e0*e0'(vI-T(u))^(-3)e0) /
%          (e0'(vI-T(u))^(-1)e0)^2
[L,d,pdflag,tmp] = toeplitzFact(-u,1);
stats(3) = stats(3) + tmp;
% finding a starting point
if pdflag % which implies u(1) < 0, can use v = 0 as a starting point
    vl = u(1);
    v = 0;
    inv_e0 = L'*(L(:,1)./d); %toeplitz(-u) \ e0;
    inv2_e0 = L'*((L*inv_e0)./d); %toeplitz(-u) \ inv_e0;
    fpv = inv2_e0(1)/inv_e0(1)-DC;
    if leq && (abs(fpv) < tol) % inequality satisfied with strict inequality
        % v = 0
        b = inv_e0/sqrt(inv_e0(1)); 
        sol = conv(b(end:-1:1),b);
        sol = sol(p+1:end);
        gphi = [u(1); 2*u(2:end)];
        return;
    end
else % use max eig-value of an augmented circulant toeplitz matrix
     % as starting point
    vl = max(0,u(1)); 
    tmp = real(fft([u;u(end-1:-1:2)])); 
    v = max(tmp);                     
    [L,d,pdflag,tmp] = toeplitzFact(-u+[v;zeros(p,1)],1);
    stats(3) = stats(3) + tmp;
    if ~pdflag, % the case when upperbound is tight, i.e. 
                % max eigs of toeplitz(u) and the augmented matrix coincide
        disp('  Wrong upperbound!!'); 
        v = v + tol; % (not sure what to use)
        %v = v + abs(u(1));
        [L,d,pdflag,tmp] = toeplitzFact(-u+[v;zeros(p,1)],1);
        stats(3) = stats(3) + tmp;
    end
    inv_e0 = L'*(L(:,1)./d); %toeplitz(-u+[v;zeros(p,1)]) \ e0;
    inv2_e0 = L'*((L*inv_e0)./d); %toeplitz(-u+[v;zeros(p,1)]) \ inv_e0;
    fpv = inv2_e0(1)/inv_e0(1) - DC;
end

if abs(fpv) < tol
    maxiter = 0; % skip iterations
else
    fv = -log(inv_e0(1))-DC*v;
    alpha = 0.01;
end

for ii = 1:maxiter % Safeguarded Newton's with backtracking
    if DEBUG, disp(ii); disp([vl v]); end
    stats(2) = stats(2) + 1;
    if ~pdflag, disp('  BUG!');  disp(ii); end
    fppv = ((-2)*inv_e0(1)*(inv_e0'*inv2_e0)+(inv2_e0(1))^2)/(inv_e0(1)^2);
    dv = -fpv/fppv;
    % line search
    t = 1;        
    while (v+t*dv) <= vl
        if DEBUG, disp('  line search bounded by lowerbound vl'); end
        t = t/2;
    end
    [L,d,pdflag,tmp] = toeplitzFact(-u+[v+t*dv;zeros(p,1)],1);
    stats(3) = stats(3) + tmp;
    while ~pdflag
        if DEBUG, disp('  line search, lowerbound vl updated'); end
        vl = v + t*dv;
        t = t/2;
        [L,d,pdflag,tmp] = toeplitzFact(-u+[v+t*dv;zeros(p,1)],1);
        stats(3) = stats(3) + tmp;
    end

    inv_e0 = L'*(L(:,1)./d); %toeplitz(-u+[v+t*dv;zeros(p,1)]) \ e0;
    while fpv < 0 && ((-log(inv_e0(1))-DC*v-fv)-t*dv*(DC + alpha*fpv) < 0)
         %((-log(inv_e0(1))-DC*(v+t*dv)) < (fv + alpha*t*fpv*dv))
        if DEBUG, disp('  line search'); end
        t = t/2;
        [L,d,pdflag,tmp] = toeplitzFact(-u+[v+t*dv;zeros(p,1)],1);
        stats(3) = stats(3) + tmp;
        inv_e0 = L'*(L(:,1)./d); %toeplitz(-u+[v+t*dv;zeros(p,1)]) \ e0;
    end
    v = v + t*dv;
    fv = -log(inv_e0(1))-DC*v;       

    inv2_e0 = L'*((L*inv_e0)./d); %toeplitz(-u+[v;zeros(p,1)]) \ inv_e0;
    fpv = inv2_e0(1)/inv_e0(1) - DC;
    if abs(fpv) < tol
        break; 
    elseif fpv > 0
        vl = v;
    end
end

if (ii <= maxiter)
else % ii == [], abs(fv) >= tol    
    ii = maxiter;
    if maxiter ~= 0
        disp('   Maximum iteration reached...');
    end
end
if leq && (v < 0), disp('   Missed case!!!'); end
if DEBUG, disp([ii v]); end
b = inv_e0/sqrt(inv_e0(1)); 
phi = -2*log(b(1));
sol = conv(b(end:-1:1),b);
sol = sol(p+1:end);
gphi = [u(1)-v; 2*u(2:end)];
stats(1) = ii; %stats(2) + stats(3);
stats(4) = cputime - stats(4);
end