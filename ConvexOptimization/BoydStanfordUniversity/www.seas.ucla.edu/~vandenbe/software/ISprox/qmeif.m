function [mu,nu,Q,W] = qmeif(U,V,Phi,Psi,tol)
% QMEIF: Matrix factorization for quadratic matrix equality and inequality.
%    [mu,nu,Q,W] = QMEIF(U,V,Phi,Psi) assumes the inputs satisfy
%       Phi(1,1)*U*U' + Phi(2,1)*U*V' + Phi(1,2)*V*U' + Phi(2,2)*V*V' = 0
%       Psi(1,1)*U*U' + Psi(2,1)*U*V' + Psi(1,2)*V*U' + Psi(2,2)*V*V' <= 0,
%    and produces a factorization for U and V such that
%       U*Q = W*diag(mu), V*Q = W*diag(nu)
%    with Q unitary and (mu(i),nu(i)) != (0,0) satisfying
%       [mu(i); nu(i)]'*Phi*[mu(i); nu(i)] = 0
%       [mu(i); nu(i)]'*Psi*[mu(i); nu(i)] <= 0, i = 1,...,r.
%    If the factorization is not successful, then (mu(i),nu(i)) == (0,0) 
%    for all i and a corresponding error message will appear.
%
%    [mu,nu,Q,W] = QMEIF(U,V,Phi,Psi,tol) sets the tolerance for the 
%    inaccuracy of the quadratic matrix equality and inequality. The 
%    default tol = 1e-10.
%
%    [mu,nu,Q,W] = QMEIF(U,V) or [mu,nu,Q,W] = QMEIF(U,V,Phi) or 
%    [mu,nu,Q,W] = QMEIF(U,V,Phi,0,tol) calls QMEF(U,V,Phi,tol), where the 
%    default values tol = 1e-10 and Phi = [1 0; 0 -1] are used if not 
%    specified in the input arguments.
%
% INPUT:
%    U, V: p x r (complex)
%    Phi, Psi: 2 x 2 Hermitian, det(Phi) < 0
%              (default Phi = [1 0; 0 -1], Psi = 0)
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
DEBUG = 0; % use 1 to print out T#, 2 to print out certain values

if nargin < 5
    tol = 1e-10; % tolerance for matrix equality and inequality
end

% no inequality Psi
if nargin < 4 || (size(Psi,1) == 1 && size(Psi,2) == 1 && Psi == 0)
    if nargin < 3
        Phi = [1 0; 0 -1];
    end
    [mu,nu,Q,W] = qmef(U,V,Phi,tol);
    return;
end

p = size(U,1);
r = size(U,2);
mu = zeros(r,1); nu = zeros(r,1); Q = zeros(r,r); W = zeros(p,r);

if size(V,1) ~= p || size(V,2) ~= r
    fprintf('Matrices U and V sizes mismatch.\n');
    return;
end
if norm(U,'fro') == 0 && norm(V,'fro') == 0
    fprintf('Zero matrices U and V.\n');
    return;
end
if size(Phi,1) ~= 2 || size(Phi,2) ~= 2 ...
        || abs(Phi(1,2)-conj(Phi(2,1))) ~= 0 ... %> tol2 ... % or ~= 0 
        || imag(Phi(1,1)) ~= 0 || imag(Phi(2,2)) ~= 0 ...
        || det(Phi) >= 0
    fprintf('Invalid Phi.\n');
    return;
end
if size(Psi,1) ~= 2 || size(Psi,2) ~= 2 ...
        || abs(Psi(1,2)-conj(Psi(2,1))) ~= 0 ... %> tol2 ... % or ~= 0
        || imag(Psi(1,1)) ~= 0 || imag(Psi(2,2)) ~= 0 % Psi ~= 0 always holds here
    fprintf('Invalid dimensions of Psi or not Hermitian.\n');
    return;
end
if norm(Phi(1,1)*(U*U') + Phi(2,1)*U*V' ...
        + Phi(1,2)*V*U' + Phi(2,2)*(V*V')) > tol
    fprintf('Quadratic matrix equation does not hold.\n');
    return;
end
if max(eig(Psi(1,1)*(U*U') + Psi(2,1)*U*V' ...
        + Psi(1,2)*V*U' + Psi(2,2)*(V*V'))) > tol 
    fprintf('Quadratic matrix inequality does not hold.\n');
    return;
end

% find R1
if norm(Phi - [0 1; 1 0],'fro') == 0
    R = eye(2);
    invRT = R;
else % other Phi = R'*[0 1; 1 0]*R
    [Qp, D] = schur(Phi);
    if DEBUG == 2, Qp, D, end
    Qp = [Qp(:,2) Qp(:,1)];
    D = abs(diag(D)); D = [D(2); D(1)];
    R = [1 1; 1 -1]*diag(sqrt(D./2)) * Qp';
    invRT = conj(R') \ eye(2);
    Psi = conj(invRT) * Psi * (conj(invRT))';
    Psi = .5*(Psi + Psi');
end

% find R
if imag(Psi(1,2)) == 0
    if DEBUG, disp('T4'); end
    if Psi(1,1) < Psi(2,2)
        % exchange the two elements of Psi and accordingly, R and invRT
        temp = Psi(2,2);
        Psi(2,2) = Psi(1,1);
        Psi(1,1) = temp;
        R = [R(2,:); R(1,:)];
        invRT = [invRT(2,:); invRT(1,:)];
        if DEBUG, disp('T5'); end
    end
else % case with complex Psi
    % Iwasaki-Hara transformation
    [Up,D] = ...
        schur([Psi(1,1) 1i*imag(Psi(1,2)); ...
               1i*imag(Psi(2,1)) Psi(2,2)]);
    if DEBUG == 2, Up, D, end
    if D(1,1) < D(2,2) 
        D = [D(2,2) 0; 0 D(1,1)]; 
        Up = [Up(:,2) Up(:,1)]; 
        if DEBUG, disp('T6'); end
    end
    % needs checking
    Up(:,1) = Up(:,1) * (abs(Up(1,1))/Up(1,1));
    Up(:,2) = Up(:,2) * (abs(Up(2,2))/Up(2,2));
    R = Up' * R;
    invRT = conj(R') \ eye(2);
    Psi = D + real(Psi(1,2))*[0 1; 1 0]; % should be Hermitian
    if DEBUG, disp('T7'); end
end
    
S = R(1,1)*U + R(1,2)*V;
T = R(2,1)*U + R(2,2)*V;
U = S;
V = T;

if Psi(2,2) > tol2 % program should not come here
    disp('Infeasible matrix inequality. Program should not come here!');
    return;
elseif Psi(1,1) <= tol2 % Psi(2,2) <= Psi(1,1) <= 0
    % redundant inequality
    [s,t,Q,W] = qmef(U,V,[0 1; 1 0],tol);
    if DEBUG, disp('T8'); end
elseif Psi(1,1) > 0 && abs(Psi(2,2)) <= tol2
    s = zeros(r,1); t = ones(r,1); W = V; Q = eye(r);
    if DEBUG, disp('T9'); end
else % Psi(2,2) < 0 < Psi(1,1)
    if DEBUG, disp('T0'); end
    % now U*V' + V*U' = 0, Psi(1,1)*U*U' <= -Psi(2,2)*V*V'
    W = V;
    V = sqrt(-Psi(2,2)/Psi(1,1)) * V;
    % now U*U' <= V*V'
    % find Lda s.t. U = V * Lda, Lda + Lda' = 0, Lda*Lda' <= I
    tmp = (V*V')-(U*U');
    if DEBUG == 2, eig(tmp)', end
    if max(eig(tmp))> tol2
        [Qu,Du] = schur(tmp,'complex');
        U2 = Qu*diag(sqrt(diag(Du)));
        U = [U U2]; V = [V zeros(size(U2))];
        if DEBUG, disp('T3'); end
    end
    % now U*V' + V*U' = 0, U*U' = V*V'
    if DEBUG == 2, disp(norm(U*V'+V*U','fro')), disp(norm(U*U'-V*V')), end
    [P,Sig,Q1] = svd(U,0);
    mm = size(U,2);
    % find Q2
    rk = rank(Sig);
    Q2 = (Sig(1:rk,1:rk) \ (P(:,1:rk)'*V))';
    if rk < mm 
        [Q22,tmp] = schur(eye(mm)-Q2*Q2','complex');
        Q2 = [Q2(:,1:rk) Q22(:,rk+1:end)];
        if DEBUG, disp('T01'); end
    end

    Lda = Q1'*Q2;
    Lda11 = Lda(1:rk,1:rk);
    if DEBUG == 2, eig(Lda11).', end
    Lda12 = Lda(1:rk,rk+1:end);
    [QL,D] = schur(Lda11,'complex');
    D = diag(D);
    Id1 = (abs(D) < 1);
    QL1 = QL(:,Id1);
    D1 = D(Id1);
    OM = diag(1 ./ sqrt(1 - D1 .* conj(D1)))*QL1'*Lda12;
    Lda = Q2 * [Lda11 Lda12; -Lda12' OM'*diag(D1')*OM] * Q2';
    
    Lda = Lda(1:r,1:r);
    [Q,rho] = schur(Lda,'complex');
    if DEBUG == 2, diag(rho), end
    W = W*Q;
    s = sqrt(-Psi(2,2)/Psi(1,1)) * diag(rho);
    t = ones(r,1);
end

uvH = (invRT)' * [s'; t'];
mu = uvH(1,:)';
nu = uvH(2,:)';
end
