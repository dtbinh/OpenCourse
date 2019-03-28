function [A,b,apx_err]=battlse(r)

% BATTLSE finds the optimal solution of the best r-term piecewise linear
% convex 'lower' approximation problem defined in the robust GP paper, 
% for given degree of approximation r.
%
% User's manual on using BATTLSE is available at: 
% http://www.stanford.edu/~boyd/rgp.html
%
% NOTE: 1. BATTLSE is designed for r<=100. 
%       2. We round each entry of A and b up to 1e-8.
%
% Author: Kan-Lin Hsiung, EE department, Stanford, October 2003.

% Subroutine called: INVAPX.

A=[]; b=[]; apx_err=[];

if (r>100)
   error('The input argument (degree of approximation) is required to be <= 100.')
end

tol=1e-14; 

k=0; 
uu=log(2); ll=0;
if r==2
   A=[0,1;1,0]; b=[0;0]; 
   uu=log(2);
else  
   while uu-ll>tol, % Find the optimal solution by bisection method.
      mid=(uu+ll)/2;
      [A,b]=invapx(mid); 
      if size(A,1)>r, ll=mid; else uu=mid; end
   end
   [A,b]=invapx(uu);
end

apx_err=uu; 

return

%-------------------------------------------------------------------------
function [F,g]=invapx(epsilon)

% INVAPX finds the optimal solution of the 'inverse' best PWL convex lower 
% approximation problem defined in the robust GP paper.
%
% INPUT: epsilon is a given error of lower approximation.
% OUTPUT: (F,g) describes the optimal solution.
%
% NOTE: epsilon < log(2) is assumed, since the cases with epsilon >= log(2) 
%       are trivial.

% Subroutine called: ADDSHP.

F=[1,0]; g=0; xys=[0,-inf]; vert=[0,log(exp(epsilon)-1)]; 
% Initialize the algorithm.

tol=1e-12; leftest_x=log(exp(epsilon)-1); M=abs(1e2*leftest_x); 
% Set some tolerant and boundary values.

while 1
   ll=-M; uu=vert(1,1); vx=vert(1,1); vy=vert(1,2); 
   while uu-ll>tol, % Find the next tangent line by bisection method.
      midx=(ll+uu)/2; midy=log(1-exp(midx)); 
      [Fi,gi]=addshp(midx,midy);
      if Fi*[vx;vy]+gi>=0, uu=midx; end
      if Fi*[vx;vy]+gi<=0, ll=midx; end
   end
   xys=[midx,midy;xys]; F=[Fi;F]; g=[gi;g]; 
   
   ll=-M; uu=xys(1,1);
   while uu-ll>tol, % Find the next vertex by bisection method.
      mid=(ll+uu)/2; midval=log(exp(mid)+exp((-gi-Fi(1)*mid)/Fi(2)));
      if midval<=epsilon, uu=mid; end
      if midval>=epsilon, ll=mid; end
   end      
   vx=ll; vy=(-gi-Fi(1)*vx)/Fi(2); 
   
   if vx>leftest_x, % Check if the stopping rule is satisfied.
      vert=[vx,vy;vert];
   else
      vert=[-gi/Fi(1),0;vert];
      F=[0,1;F]; g=[0;g]; xys=[-inf,0;xys];
      break
   end   
end

S=1e8; F=round(F*S)/S; g=round(g*S)/S; 
% Rounding each entry of F and g up to 1e-8.

return

%-------------------------------------------------------------------------
function [Fi,gi]=addshp(x,y)

% Find the tangent line of lse(x,y)=0 at given point (x,y). 
fgrad=[exp(x);exp(y)]/sum(exp(x)+exp(y)); 
Fi=fgrad'; gi=-fgrad'*[x;y];
