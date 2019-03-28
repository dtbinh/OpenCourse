function [x1,x2,x3,x4,v] = decoding

% data for exercise 1.19 (Fall 2015)
% decoding using inner products

randn('state',0);
t = linspace(0,2*pi,200)';
X = [2*sin(5*t).*exp(-t/2) cos(4*t) sin(2*t)+1  -cos(4*t)];
v = X(:,2)/3 + 0.05*randn(200,1);
x1 = X(:,1);  x2 = X(:,2);  x3 = X(:,3);  x4 = X(:,4);

