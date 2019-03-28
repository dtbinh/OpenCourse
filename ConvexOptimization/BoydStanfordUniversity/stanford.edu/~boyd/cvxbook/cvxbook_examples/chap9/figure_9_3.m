% Generates figures 9.3, 9.4, 9.5, 9.11, 9.12, 9.13, 9.14, 9.15, 
% 9.19, 9.20 of Boyd & Vandenberghe, Convex Optimization
%
% Gradient method, steepest descent with quadratic norms, and Newton 
% method for an example with two variables.
%

% f(x) = sum(exp(A*x+b)) 
A = [1 3; 1 -3; -1 0];
b = -0.1*[1; 1; 1];

% backtracking parameters
alpha = 0.1;  
beta = 0.7;

% starting point
x0 = [-1;1];
f0 = sum(exp(A*x0+b));

% Calculate the minimum.
x = x0;
for i=1:10
    grad = A' * exp(A*x+b);
    hess = A' * diag(exp(A*x+b)) * A;
    x = x - hess\grad;
end;
xmin = x;
fmin = sum(exp(A*xmin+b));


% Gradient descent with backtracking line search. 

figure(1) 

% Figure 9.3: show five contour lines: f0+k*Delta, k=-3,-2,-1,0,1 where 
% Delta = (f0-fmin)/4.
nangles = 1000;
thetas = linspace(0,2*pi,nangles);
cont1 = zeros(2,nangles);
cont2 = zeros(2,nangles);
cont3 = zeros(2,nangles);
cont4 = zeros(2,nangles);
Delta = (f0-fmin)/4;
for i=1:nangles
    nopts = 1000;
    ray  = [ linspace(0,4,nopts)*cos(thetas(i)); 
        linspace(0,4,nopts)*sin(thetas(i)) ];
    fray = sum(exp(A*ray+b(:,ones(1,nopts))));
    ind = max(find(fray < f0));
    cont1(:,i) = (ray(:,ind) + ray(:,ind+1))/2;
    ind = max(find(fray < f0 - 3*Delta));
    cont2(:,i) = (ray(:,ind) + ray(:,ind+1))/2;
    ind = max(find(fray < f0 - 2*Delta));
    cont3(:,i) = (ray(:,ind) + ray(:,ind+1))/2;
    ind = max(find(fray < f0 - Delta));
    cont4(:,i) = (ray(:,ind) + ray(:,ind+1))/2;
    ind = max(find(fray < f0 + Delta));
    cont5(:,i) = (ray(:,ind) + ray(:,ind+1))/2;
end;
plot(cont1(1,:), cont1(2,:), '--', cont2(1,:), cont2(2,:), '--', ...
    cont3(1,:), cont3(2,:), '--', cont4(1,:), cont4(2,:), '--', ...
    cont5(1,:), cont5(2,:), '--');
axis('equal')
axis off
hold on

% Run the gradient method.
maxiters = 100;
xs = zeros(2,maxiters);
fs = zeros(1,maxiters);
x = x0;
for i=1:maxiters
    f = sum(exp(A*x+b));
    xs(:,i) = x;
    fs(i) = f;
    if (f-fmin < 1e-10),  break; end;
    g = A'*exp(A*x+b);
    v = -g;
    s = 1;
    for k=1:10 
        xnew = x+s*v;
        fxnew = sum(exp(A*xnew+b));
        if (fxnew < f + s*alpha*g'*v) break;  else s = s*beta; end;
    end;
    x = x + s*v;
end;

% Finish figure 9.3.
plot(xs(1,1:i), xs(2,1:i), 'o'); 
plot(xs(1,1:i), xs(2,1:i), '-'); 
hold off

% Figure 9.4: error vs. iteration for backtracking line search.
figure(2)
semilogy([0:i-1], fs(1:i)-fmin, '-', [0:i-1], fs(1:i)-fmin, 'o');
j=16;  text(j-1, fs(j)-fmin, 'backtracking');
axis([0 i-1 1e-15 1e5]);


% Gradient descent with exact line search.

% Figure 9.5: show the five contourlines.
figure(3) 
plot(cont1(1,:), cont1(2,:), '--', cont2(1,:), cont2(2,:), '--',...
     cont3(1,:), cont3(2,:), '--', cont4(1,:), cont4(2,:), '--',...
     cont5(1,:), cont5(2,:), '--');
axis('equal')
axis off
hold on

% Run the gradient method, using bisection in the exact line search.
maxiters = 100;
xs = zeros(2,maxiters);
fs = zeros(1,maxiters);
x = x0;
for i=1:maxiters
    f = sum(exp(A*x+b));
    xs(:,i) = x;
    fs(i) = f;
    g = A'*exp(A*x+b);
    if (f-fmin < 1e-10), break; end;
    v = -g;
    % exact line search by bisection
    l = 0; u = 1;
    gl = v'*g;  
    gr = v'*A'*exp(A*(x+u*v)+b);
    while (gr<0), 
        u = 2*u; 
        gr = v'*A'*exp(A*(x+u*v)+b);
    end;
    while (u-l > 1e-3)
        s = (u+l)/2;
        gs = v'*A'*exp(A*(x+s*v)+b);
        if (gs<0), l=s; 
        elseif (gs>0), u=s; 
        else, break, end;
    end;
    x = x + s*v;
end;

% Finish figure 9.5.
plot(xs(1,1:i), xs(2,1:i), 'o'); 
plot(xs(1,1:i), xs(2,1:i), '-'); 
hold off

% Finish figure 9.4.
figure(2) 
hold on 
semilogy([0:i-1], fs(1:i)-fmin, '-', [0:i-1], fs(1:i)-fmin, 'o');
j=16;  text(j-1, fs(j)-fmin, 'exact line search');
xlabel('k');  
ylabel('error');
hold off


% Steepest descent in first quadratic norm.
H = [2 0; 0 8];

% Figure 9.11: show the five contour lines.
figure(4)
plot(cont1(1,:), cont1(2,:), '--', cont2(1,:), cont2(2,:), '--', ...
     cont3(1,:), cont3(2,:), '--', cont4(1,:), cont4(2,:), '--', ...
     cont5(1,:), cont5(2,:), '--');
axis('equal')
axis off
hold on

% Run the steepest descent method.
maxiters = 100;
xs = zeros(2,maxiters);
fs = zeros(1,maxiters);
x = x0;
for i=1:maxiters
    f = sum(exp(A*x+b));
    xs(:,i) = x;
    fs(i) = f;
    if (f-fmin < 1e-10) break; end;
    g = A'*exp(A*x+b);
    v = -H\g;
    s = 1;
    for k=1:10 
        xnew = x+s*v;
        fxnew = sum(exp(A*xnew+b));
        if (fxnew < f + s*alpha*g'*v), break;  else, s = s*beta; end;
    end;
    x = x + s*v;
end;
plot(xs(1,1:i), xs(2,1:i), 'o'); 
plot(xs(1,1:i), xs(2,1:i), '-'); 

% Plot ellipses around x0 and x1.
x = x0;
ell = x(:,ones(1,nangles)) + sqrtm(H)\[cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
x = xs(:,2);
ell = x(:,ones(1,nangles)) + sqrtm(H)\[cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
hold off

% Figure 9.14: same plot after change of coordinates.
figure(5)
cont1sc = sqrtm(H)*cont1;
cont2sc = sqrtm(H)*cont2;
cont3sc = sqrtm(H)*cont3;
cont4sc = sqrtm(H)*cont4;
cont5sc = sqrtm(H)*cont5;
plot(cont1sc(1,:), cont1sc(2,:), '--', cont2sc(1,:), cont2sc(2,:), ...
    '--', cont3sc(1,:), cont3sc(2,:), '--', cont4sc(1,:), ...
    cont4sc(2,:), '--', cont5sc(1,:), cont5sc(2,:), '--');
axis('equal')
axis off
hold on
xssc = sqrtm(H)*xs;
plot(xssc(1,1:i), xssc(2,1:i), 'o'); 
plot(xssc(1,1:i), xssc(2,1:i), '-'); 

% Plot ellipses around x0 and x1.
x = x0;
ell = sqrtm(H) * x(:,ones(1,nangles)) + [cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
x = xs(:,2);
ell = sqrtm(H) * x(:,ones(1,nangles)) + [cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
hold off

% Figure 9.13: error vs. iteration for steepest descent in the first
% norm.
figure(6)
semilogy([0:i-1], fs(1:i)-fmin, '-', [0:i-1], fs(1:i)-fmin, 'o');
text(10, fs(11)-fmin,'P1');
hold on

% Steepest descent in second quadratic norm.
H = [8 0; 0 2];

% Figure 9.12: show the five contour lines.
figure(7)
plot(cont1(1,:), cont1(2,:), '--', cont2(1,:), cont2(2,:), '--', ...
    cont3(1,:), cont3(2,:), '--', cont4(1,:), cont4(2,:), '--', ...
    cont5(1,:), cont5(2,:), '--');
axis('equal')
axis off
hold on

% Run the steepest descent method.
maxiters = 100;
xs = zeros(2,maxiters);
fs = zeros(1,maxiters);
x = x0;
for i=1:maxiters
    f = sum(exp(A*x+b));
    xs(:,i) = x;
    fs(i) = f;
    if (f-fmin < 1e-10) break; end;
    g = A'*exp(A*x+b);
    v = -H\g;
    s = 1;
    for k=1:10 
        xnew = x+s*v;
        fxnew = sum(exp(A*xnew+b));
        if (fxnew < f + s*alpha*g'*v) break;  else s = s*beta; end;
    end;
    x = x + s*v;
end;
plot(xs(1,1:i), xs(2,1:i), 'o'); 
plot(xs(1,1:i), xs(2,1:i), '-'); 

% Plot ellipses around x0 and x1.
x = x0;
ell = x(:,ones(1,nangles)) + sqrtm(H) \ [cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
x = xs(:,2);
ell = x(:,ones(1,nangles)) + sqrtm(H) \ [cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
hold off

% Figure 9.15: same plot after change of coordinates.
figure(8)
cont1sc = sqrtm(H)*cont1;
cont2sc = sqrtm(H)*cont2;
cont3sc = sqrtm(H)*cont3;
cont4sc = sqrtm(H)*cont4;
cont5sc = sqrtm(H)*cont5;
plot(cont1sc(1,:), cont1sc(2,:), '--', cont2sc(1,:), cont2sc(2,:), ...
    '--', cont3sc(1,:), cont3sc(2,:), '--', cont4sc(1,:), ...
    cont4sc(2,:), '--', cont5sc(1,:), cont5sc(2,:), '--');
axis('equal')
axis off
hold on
xssc = sqrtm(H)*xs;
plot(xssc(1,1:maxiters), xssc(2,1:maxiters), 'o'); 
plot(xssc(1,1:maxiters), xssc(2,1:maxiters), '-'); 

% Plot ellipses around x0 and x1.
x = x0;
ell = sqrtm(H) * x(:,ones(1,nangles)) + [cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
x = xs(:,2);
ell = sqrtm(H) * x(:,ones(1,nangles)) + [cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
hold off


% Finish figure 9.13: error vs. iteration for steepest descent in the 
% second  norm.

figure(6)
semilogy([0:i-1], fs(1:i)-fmin, '-', [0:i-1], fs(1:i)-fmin, 'o');
text(30, fs(31)-fmin,'P2');
axis([0 40 1e-15 1e5]);
xlabel('k');  ylabel('error');
hold off


% Newton's method

% Figure 9.19.  Plot the five contour lines.
figure(9) 
plot(cont1(1,:), cont1(2,:), '--', cont2(1,:), cont2(2,:), '--', ...
    cont3(1,:), cont3(2,:), '--', cont4(1,:), cont4(2,:), '--', ...
    cont5(1,:), cont5(2,:), '--');
axis('equal')
axis off
hold on

% Run Newton's method.
maxiters = 10;
xs = zeros(2,maxiters);
fs = zeros(1,maxiters);
steps = zeros(1,maxiters-1);
x = x0;
for i=1:maxiters
    f = sum(exp(A*x+b));
    xs(:,i) = x;
    fs(i) = f;
    g = A'*exp(A*x+b);
    H = A'*diag(exp(A*x+b))*A;
    v = -H\g;
    if (abs(v'*g) < 1e-8), break; end;
    s = 1;
    for k=1:10 
        xnew = x+s*v;
        fxnew = sum(exp(A*xnew+b));
        if (fxnew < f + s*alpha*g'*v), break;  else, s = s*beta; end;
    end;
    steps(i) = s;
    x = x + s*v;
end;
plot(xs(1,1:i), xs(2,1:i), 'o'); 
plot(xs(1,1:i), xs(2,1:i), '-'); 

% Plot ellipses around x0 and x1.
x = x0;
H = A'*diag(exp(A*x+b))*A;
ell = x(:,ones(1,nangles)) + sqrtm(H) \ [cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
x = xs(:,2);
H = A'*diag(exp(A*x+b))*A;
ell = x(:,ones(1,nangles)) + sqrtm(H) \  [cos(thetas); sin(thetas)];
plot(ell(1,:), ell(2,:), '-');
hold off

% Figure 9.20.  Error vs. iteration.
figure(10)
semilogy([0:i-1], fs(1:i)-fmin, '-', [0:i-1], fs(1:i)-fmin, 'o');
axis([0 i-1 1e-15 1e5]);
xlabel('k');
ylabel('error');
