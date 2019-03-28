function [mu,nu,Q,W] = qmef(U,V,Phi,tol)
% QMEF: Matrix factorization for quadratic matrix equality.
%    [mu,nu,Q,W] = QMEF(U,V,Phi) assumes the inputs satisfy
%       Phi(1,1)*U*U' + Phi(2,1)*U*V' + Phi(1,2)*V*U' + Phi(2,2)*V*V' = 0,
%    and produces a factorization for U and V such that
%       U*Q = W*diag(mu), V*Q = W*diag(nu)
%    with Q unitary and (mu(i),nu(i)) != (0,0) satisfying
%       [mu(i); nu(i)]'*Phi*[mu(i); nu(i)] = 0, i = 1,...,r.
%    If the factorization is not successful, then (mu(i),nu(i)) == (0,0) 
%    for all i and a corresponding error message will appear.
%
%    [mu,nu,Q,W] = QMEF(U,V) assumes and uses Phi = [1 0; 0 -1].
%
%    [mu,nu,Q,W] = QMEF(U,V,Phi,tol) sets the tolerance for the inaccuracy
%    of the quadratic matrix equation. The default tol = 1e-10.
%
% INPUT: 
%    U, V: p x r (complex)
%    Phi: 2 x 2 Hermitian , det(Phi) < 0  (default Phi = [1 0; 0 -1])
%    tol : tolerance of inaccuracy, real positive scalar  
%          (default tol = 1e-10)
% OUTPUT:
%    mu, nu: r x 1 (complex)
%    Q: r x r unitary
%    W: p x r (complex)

% Procedure follows lemma 1 and lemma 2 in the paper:
%    Semidefinite representations of gauge functions for structured 
%    low-rank matrix decomposition
%    (Available at https://arxiv.org/abs/1604.02500)
% Author: Hsiao-Han (Suzi) Chao (suzihchao@ucla.edu)
% Date: 11/15/2016

% The following two parameters can be modified for debugging
tol2 = 1e-10; % numerical accuracy
DEBUG = 0; % use 1 to print out T#, use 2 to print out certain values

% tolerance for matrix equality
if nargin < 4
    tol = 1e-10; 
end
p = size(U,1);
r = size(U,2);
mu = zeros(r,1); nu = zeros(r,1); Q = zeros(r,r); W = zeros(p,r);

if size(V,1) ~= p || size(V,2) ~= r
    fprintf('Matrices U and V sizes mismatch.\n');
    return;
end

if nargin < 3
    % Phi == [1 0; 0 -1];
    R = eye(2);
else % other Phi
    if size(Phi,1) ~= 2 || size(Phi,2) ~= 2 ...
            || abs(Phi(1,2)-conj(Phi(2,1))) > tol2 ... % or ~= 0
            || imag(Phi(1,1)) ~= 0 || imag(Phi(2,2)) ~= 0 ...
            || det(Phi) >= 0
        fprintf('Invalid Phi.\n');
        return;
    end
    
    [Qp, D] = schur(Phi,'complex');
    Qp = [Qp(:,2) Qp(:,1)]; 
    D = abs(diag(D)); D = [D(2); D(1)];
    R = diag(sqrt(D)) * Qp';
    S = R(1,1)*U + R(1,2)*V;
    T = R(2,1)*U + R(2,2)*V;
    U = S;
    V = T;
end

if norm(U*U'-V*V') > tol
    if DEBUG, fprintf('norm(UU''-VV'') = %f\n',norm(U*U'-V*V')); end
    fprintf('Quadratic matrix equation does not hold.\n');
    return;
end

% compute Q1
[P,Sig,Q1] = svd(U,0);
% find Q2
rk = rank(Sig);
Q2 = (Sig(1:rk,1:rk) \ (P(:,1:rk)'*V))';
if DEBUG, disp('T1'); end
if rk < r
    [Q22,tmp] = schur(eye(r)-Q2*Q2','complex');
    Q2 = [Q2(:,1:rk) Q22(:,rk+1:end)];
    if DEBUG, disp('T2'); end
end

[Q, S] = schur(Q2*Q1','complex');
s = diag(S);
W = V*Q;
invR_ = conj(R) \ eye(2);
uvH = invR_ * [s'; ones(1,length(s))];
mu = uvH(1,:)';
nu = uvH(2,:)';
end

