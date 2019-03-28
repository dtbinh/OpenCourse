function [A, Ides] = illumdata 

% illumination problem

% lamp positions 
L = [0.1 0.3 0.4 0.6 0.8 0.9 0.95
     1.0 1.1 0.6 0.9 0.9 1.2 1.00];
m = size(L,2);    % number of lamps

% begin and endpoints of patches 
V = [linspace(0,1,12);
     0   0.1 0.2 0.2 0.1 0.2 0.3 0.2 0   0   0.2 0.1];
n = size(V,2)-1;  % number of patches

% desired illumination
Ides = 2;

% construct A
dV = V(:,2:n+1)-V(:,1:n);    % tangent to patches
VI = V(:,1:n)+.5*dV;         % midpoint of patches
A = zeros(n,m);
for i=1:n
  for j=1:m
    dVI = L(:,j)-VI(:,i);  
    dVperp = null(dV(:,i)');  % upward pointing normal 
    if dVperp(2)<0
      dVperp = -dVperp;
    end
    A(i,j) = max(0,dVI'*dVperp/(norm(dVI)*norm(dVperp)))./norm(dVI)^2;
  end
end
