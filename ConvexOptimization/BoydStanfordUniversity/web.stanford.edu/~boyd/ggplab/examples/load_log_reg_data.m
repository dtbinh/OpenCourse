% Data for the GP logistic regressioni example.
% Obtained from the book Boyd & Vandenberghe "Convex Optimization"
% (figure 7.1, page 355)
randn('state',0);
rand('state',0);

a =  1;
b = -5;

m = 100;
u = 10*rand(m,1);
y = (rand(m,1) < exp(a*u+b)./(1+exp(a*u+b)));
