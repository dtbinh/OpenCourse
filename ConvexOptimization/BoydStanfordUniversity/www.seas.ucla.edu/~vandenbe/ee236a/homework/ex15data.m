function [A, b] = ex15data(s) 

% [A, b] = ex15data('small') or  [A, b] = ex15data('large')
%
% Exercise 15, Fall 2013 (placement).
 
% N = number of connections
% n = number of cells
% m = number of I/O terminals

if (strcmp(s,'small'))
   N=150; n=50; m=36; sd=3;           
elseif (strcmp(s,'large'))
   N=300; n=100; m=66; sd = 602;     
else
   disp('Invalid input argument.');  
   return;
end;

rand('seed',sd);
randn('seed',sd);

A = zeros(N,n);
b = zeros(N,1);

for i=1:N-m
   ind1 = 1+round((n-1)*rand);
   ind2 = 1+round((n-1)*rand);
   while (ind1==ind2), ind2 = 1+round((n-1)*rand); end;
   A(i,min(ind1,ind2)) = 1;
   A(i,max(ind1,ind2)) = -1;
end;
for i=N-m+[1:m/2]
   ind1 = 1+round((n-1)*rand);
   A(i,ind1) = 1;
   b(i)=-1;
end;
for i=N-m/2+[1:m/2]
   ind1 = 1+round((n-1)*rand);
   A(i,ind1) = 1;
   b(i)=1;
end;
