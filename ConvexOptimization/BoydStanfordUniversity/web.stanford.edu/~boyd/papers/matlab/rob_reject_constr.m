function rob_reject_constr()
% Matlab script for the paper "Array Signal Processing with
% Robust Rejection Constraints via Second-Order Cone Programming",
% by A. Mutapcic, S.-J. Kim, and S. Boyd.
% In Proceedings of the 40th Asilomar Conference on Signals, Systems,
% and Computers, Oct. 29 - Nov. 1, 2006, Pacific Grove, CA.
% 
% This script requires CVX software, which is available at
% http://www.stanford.edu/~boyd/cvx/

% set random seed to repeat experiments
rand('state',1);
randn('state',0);

% angle sampling 
theta  = [1:1:360]';

%********************************************************************
% antenna specs
%********************************************************************
lambda = 1;           % wavelength
theta_tar = 45;       % desired (target) direction

% interference directions
theta_rej = [10 65 120 200];
% rejection level
EPS = 0.05; % 20*log10(0.05) = -26 dB

% random array of n antenna elements
% (uniformly distributed on [0,L]-by-[0,L] square)
n = 20;
L = 5;
loc = L*rand(n,2);

%********************************************************************
% construct optimization data
%********************************************************************
% build matrix A that relates w and y(theta), ie, y = A*w
A = kron(cos(pi*theta/180), loc(:,1)') + kron(sin(pi*theta/180), loc(:,2)');
A = exp(2*pi*i/lambda*A);

% desired angle constraint
Atar = A(theta_tar,:);

% stopband constraint matrix
Areject = [];
Arej = {};
Prej = {};

for k = 1:length(theta_rej)
  Areject = [Areject; A(theta_rej(k),:)];
  Arej{k} = A(theta_rej(k),:)';
  Prej{k} = eye(n);
end


%********************************************************************
% nominal optimal solution
%********************************************************************
cvx_begin
  variable wnom(n) complex
  minimize( norm(wnom) )
  subject to
    real( Atar*wnom ) >= 1;
    for k = 1:length(Arej)
      abs( Arej{k}'*wnom ) <= EPS;
    end
cvx_end

% check if problem was successfully solved
disp(['Problem is ' cvx_status])
if ~strcmp(cvx_status,'Solved')
  error('Check problem solution...')
end

fprintf(1,'The nominal max interference angle gain is %3.2f dB.\n\n',...
          20*log10( max(abs( Areject*wnom )) ) );


wc_rej_array = [];
wc_rej_array_nom = [];
obj_value_array = [];

rho = linspace(0,0.2,20);
for r = 1:length(rho)
  fprintf(1,'processing rho = %f\n',rho(r))
  %********************************************************************
  % robust optimal solution
  %********************************************************************
  cvx_begin
    variable w(n) complex
    minimize( norm(w) )
    subject to
      real( Atar*w ) >= 1;

      for k = 1:length(Arej)
        sup_ellip_abs( w, Arej{k}, rho(r)*Prej{k} ) <= EPS;
      end

  cvx_end

  % check if problem was successfully solved
  disp(['Problem is ' cvx_status])
  if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
    disp('Check problem solution...')
    continue
  end

  obj_value_array = [obj_value_array norm(w)];

  fprintf(1,'The max interference angle gain is %3.2f dB.\n\n',...
            20*log10( max(abs( Areject*w )) ) );

  %********************************************************************
  % check the worst-case performance
  %********************************************************************
  wc_rej     = [];
  wc_rej_nom = [];

  cvx_quiet(true);
  for k = 1:length(Arej)
    wc_rej(k,1)     = sup_ellip_abs( w, Arej{k}, rho(r)*Prej{k} );
    wc_rej_nom(k,1) = sup_ellip_abs( wnom, Arej{k}, rho(r)*Prej{k} );
  end
  cvx_quiet(false);

  wc_rej_array     = [wc_rej_array max(wc_rej)];
  wc_rej_array_nom = [wc_rej_array_nom max(wc_rej_nom)];

end % for loop

%********************************************************************
% plots
%********************************************************************
% worst-case rejection level versus rho
figure(1);
plot(rho,obj_value_array,rho,norm(wnom)*ones(length(rho),1),'r--');
xlabel('rho'), ylabel('Pw')
set(gca,'FontSize',14);
axis([0 0.2 0.235 0.25]), grid on,

% objective value versus rho
figure(2);
plot(rho,wc_rej_array,rho,wc_rej_array_nom,'r--');
xlabel('rho'), ylabel('wc-null')
grid on
set(gca,'FontSize',14);

%********************************************************************
% helper function for the robust gain constraint
%********************************************************************
function cvx_optval = sup_ellip_abs(w,a_nom,P)
% SUP_ELLIP_ABS robust absolute value function over an ellipsoid.
%
% Gives an epigraph formulation of the following convex function g(w):
%
%   g(w) = sup_{ a \in E } |w*a| 
%
% where ellipsoid E = { Pu + a_nom | ||u|| <= 1 }.

cvx_optval = abs(a_nom'*w) + norm(P'*w);
