%% Basic Example

cvx_begin
  variables x y 
  minimize(x + y)
  subject to
    x >= 1
    y == 2
cvx_end

%% LP Example

A = randn(3, 10);
b = A * rand(10, 1);
c = -rand(10, 1);

cvx_begin
  variables x(10)
  maximize(c' * x)
  subject to
    A * x == b
    x >= 0
cvx_end

%% SDP Example

C = randn(3);
A = C' * C - rand(3);

cvx_begin sdp
  variable X(3, 3)
  minimize(norm(A - X))
  subject to
    X >= 0
cvx_end

%% Assignment Example

A = randn(3,3);

cvx_begin
  variable x(3)
  OBJ = 0
  for i = 1:3
    OBJ = OBJ + norm(A(:,i) - x);
  end
  minimize(OBJ)
cvx_end

