function [A,B,xy]= cut_grid_data;
% generate a cut-grid graph for the ICM 2006 talk example
% graph has 64 nodes and 95 edges
% A is an n x m incidence matrix (n is number of nodes, m is number of edges)
% B is an n x n adjacency matrix
% xy is the location data
%
% Original code by Arpita Ghosh, modified for ICM06 talk by Almir Mutapcic
n = 8; 
r1 = sqrt(2); 
x= zeros(1,n^2); 
y = zeros(1,n^2);
for i=1:n
    x(n*(i-1) + 1 : n*i) = 1:n; 
end
for j = 1:n
    y(n*(j-1)+1:n*j) = j*ones(1,n); 
end
xy = [x' y']; 
Dist1 = zeros(n^2);
for i=1:n^2;
    for j=i:n^2;
        Dist1(i,j) = norm( xy(i,:) - xy(j,:) );
        X(i,j) = 0.5*(xy(i,1)+xy(j,1)); 
        Y(i,j) = 0.5*(xy(i,2)+xy(j,2)); 
    end;
end;
Dist1 = Dist1 + Dist1';

% find the adjacency matrix
Adj1 = Dist1 < r1;
Adj1 = Adj1 - diag(diag(Adj1));
Adj1(41,49)  = 0 ; Adj1(49,41)  = 0;
Adj1(42,50) = 0; Adj1(50,42)  = 0; 
Adj1(8,16) = 0; Adj1(16,8) = 0;
Adj1(16,24) = 0; Adj1(24,16) = 0;
Adj1(7,15) = 0; Adj1(15,7) = 0;
Adj1(15,23) = 0; Adj1(23,15) = 0;
Adj1(1,10) = 0; Adj1(10,1) = 0; 
Adj1(13,21) = 0; Adj1(21,13) = 0;
Adj1(5,13) = 0; Adj1(13,5) = 0; 
Adj1(14,22) = 0; Adj1(22,14) = 0;
Adj1(6,14) = 0; Adj1(14,6) = 0; 
Adj1(43,51) = 0; Adj1(51,43)  = 0; 
Adj1(44,52) = 0; Adj1(52,44)  = 0; 
Adj1(45,53) = 0; Adj1(53,45)  = 0; 
Adj1(46,54) = 0; Adj1(54,46)  = 0; 

Adj1(41,42)  = 0 ; Adj1(42,41)  = 0;
Adj1(33,34)  = 0 ; Adj1(34,33)  = 0;
Adj1(25,26)  = 0 ; Adj1(26,25)  = 0;

X = X.*Adj1; 
Y = Y.*Adj1; 
X = X+X'; 
Y = Y+Y'; 
m = sum(sum(Adj1))/2;

Inc1 = zeros(n^2,m);
l = 0;
for i=1:(n^2-1);
    for j=i+1:n^2;
        if Adj1(i,j)>0.0001
            l = l + 1;
            Inc1(i,l) =  1;
            Inc1(j,l) = -1;
        end;
    end;
end;

B = Adj1;
A = Inc1;
