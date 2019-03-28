% Matlab script for the paper "Beamforming with Uncertain Weights",
% by A. Mutapcic, S.-J. Kim, and S. Boyd.
%
% We consider beamformers with complex multiplicative uncertainty
% in the element weights, and magnitude bound on each uncertainty.
% In this case, the robust optimal beamformers can be designed using
% the weighted complex l_1-regularization of the nominal problem.
%
% Nominal problem: minimize the rejection level over the rejection
%                  band, subject to unit array gain in the desired
%                  direction theta_des.
%
%   minimize   max |y(theta)|          for theta in rejection band 
%       s.t.   Re y(theta_des) >= 1
%
% This script requires CVX software, which is available at
% http://www.stanford.edu/~boyd/cvx/

%********************************************************************
% antenna specs for our problem instance
%********************************************************************
lambda = 1;           % wavelength
theta_des = 60;       % desired direction (should be an integer -- discretization)
half_beamwidth = 20;  % half beamwidth around the target direction

% rectangular 6 x 6 grid of 36 sensor (antenna) elements
m = 6; n = m^2;
d = 0.45*lambda;

% sensor locations
loc = zeros(n,2);
for x = 0:m-1
  for y = 0:m-1
    loc(m*y+x+1,:) = [x y];
  end
end
loc = loc*d;

% angle sampling 
theta  = [1:1:360]';

%********************************************************************
% construct optimization data
%********************************************************************
% build matrix A that relates weights w and composite array output y = A*w
A = kron(cos(pi*theta/180), loc(:,1)') + kron(sin(pi*theta/180), loc(:,2)');
A = exp(2*pi*i/lambda*A);

% desired direction constraint matrix
[diff_closest, ind_closest] = min( abs(theta - theta_des) );
Ades  = A(ind_closest,:);
thdes = theta(ind_closest);

% rejection band constraint matrix
ind = find(theta <= (theta_des-half_beamwidth) | ...
           theta >= (theta_des+half_beamwidth) );
Astop = A(ind,:);
thstop = theta(ind);

%********************************************************************
% nominal beamforming problem
%********************************************************************
cvx_begin
  variable w_nom(n) complex
  minimize( max( abs( Astop*w_nom ) ) )
  subject to
    real( Ades*w_nom ) >= 1;
cvx_end

% check if problem was successfully solved
disp(['Problem is ' cvx_status])
if ~strcmp(cvx_status,'Solved')
  return
end

% compute the optimal rejection level
mrl_nom = max(abs(Astop*w_nom));
fprintf(1,'The minimum rejection level is %3.2f dB.\n\n', 20*log10(mrl_nom) );

%********************************************************************
% robust beamforming problem 
%********************************************************************
% regularization constants (they are unity in this example)
mu = max( abs( Astop ) );
nu = abs( Ades );

% loop over different rhos and recompute robust optimal solution
rho = linspace(0,0.15,40);

mrl_rob    = [];
mrl_tik_wc = [];
mrl_nom_wc = [];

cvx_quiet(true);

for k = 1:length(rho)
  % optimal robust solution
  cvx_begin
    variable w_rob(n) complex
    minimize( max( abs( Astop*w_rob ) ) + rho(k)*mu*abs(w_rob) )
    subject to
      real( Ades*w_rob ) >= 1 + rho(k)*nu*abs(w_rob);
  cvx_end

  % check if problem was successfully solved
  disp(['rho = ' num2str(rho(k))])
  disp(['  Optimal worst-case robust problem is ' cvx_status '.'])
  if ~strcmp(cvx_status,'Solved')
    return
  end

  % worst-case robust performance
  mrl_rob(k) = 20*log10( max(abs(Astop*w_rob)) + rho(k)*mu*abs(w_rob) );

  % Tikhonov l_2-regularized robust solution
  % we have found that MU = 2 gives excellent performance (experimental)
  MU = 2; 
  cvx_begin
    variable w_tik(n) complex
    minimize( max( abs( Astop*w_tik ) ) + MU*rho(k)*sum_square_abs(w_tik) )
    subject to
      % note that the solution will not satisfy this constraint robustly
      real( Ades*w_tik ) >= 1;
  cvx_end

  % check if problem was successfully solved
  disp(['  Tikhonov l_2-regularized problem  is ' cvx_status '.'])
  if ~strcmp(cvx_status,'Solved')
    return
  end

  % scale Tikhonov solution so it satisfies real( Ades*w_tik ) == 1
  w_tik = w_tik / ( real(Ades*w_tik) - rho(k)*nu*abs(w_tik) );

  % worst-case Tikhonov performance
  mrl_tik_wc(k) = 20*log10( max(abs(Astop*w_tik)) + rho(k)*mu*abs(w_tik) );

  % also compute the worst-case performance of the nominal design
  % w_nom = w_nom / ( real(Ades*w_nom) - rho(k)*nu*abs(w_nom) );
  mrl_nom_wc(k) = 20*log10( max(abs(Astop*w_nom)) + rho(k)*mu*abs(w_nom)  );
end

cvx_quiet(false);

%********************************************************************
% plots
%********************************************************************
figure(1), clf
plot(rho,mrl_nom_wc,'r--','LineWidth',1.15), hold on,
plot(rho,mrl_tik_wc,'b--','LineWidth',1.15), hold on,
plot(rho,20*log10(mrl_nom + rho./(1-rho)),'r:','LineWidth',1.15), hold on,
plot(rho,mrl_rob,'-','Color',[0 .5 0],'LineWidth',1.15), hold off
xlabel('rho'), ylabel('mrl')
axis([0.01 .15 -25 -10])

% build matrix A for plotting
A = kron(cos(pi*theta/180), loc(:,1)') + kron(sin(pi*theta/180), loc(:,2)');
A = exp(2*pi*i/lambda*A);

% beamformer plots
figure(2), clf
Rho = rho(end);
deg = 20;
delta_rand = (2*Rho*rand(n,1)-1 + i*(2*Rho*rand(n,1)-1));
delta_nom  = Rho*exp(i*angle(Astop(deg,:)*w_nom)/n)*...
             (Astop(deg,:)'.*w_nom)./abs(Astop(deg,:)'.*w_nom);
delta_rob  = Rho*exp(i*angle(Astop(deg,:)*w_rob)/n)*...
             (Astop(deg,:)'.*w_rob)./abs(Astop(deg,:)'.*w_rob);
delta_tik  = Rho*exp(i*angle(Astop(deg,:)*w_tik)/n)*...
             (Astop(deg,:)'.*w_tik)./abs(Astop(deg,:)'.*w_tik);

y_nom = abs( A*w_nom );
y_nom_wc = abs( A*(w_nom.*(1 + delta_nom)) );
y_rob_wc = abs( A*(w_rob.*(1 + delta_rob)) );
y_tik_wc = abs( A*(w_tik.*(1 + delta_tik)) );

dBYnom = 20*log10(y_nom/y_nom(thdes));
dBYnom_wc = 20*log10(y_nom_wc/y_nom_wc(thdes));
dBYrob_wc = 20*log10(y_rob_wc/y_rob_wc(thdes));
dBYtik_wc = 20*log10(y_tik_wc/y_tik_wc(thdes));

plot(theta,dBYnom_wc,'r--','LineWidth',1.15), hold on
% plot(theta,dBYtik_wc,'b:','LineWidth',1.15), hold on
plot(theta,dBYrob_wc,'-','Color',[0 .5 0],'LineWidth',1.15), hold on

ymin = -60; ymax = 20;
plot([thdes+half_beamwidth thdes+half_beamwidth],[ymin ymax],'k:',...
     [thdes thdes],[ymin ymax],'k:',...
     [thdes-half_beamwidth thdes-half_beamwidth],[ymin ymax],'k:')
hold off
axis([0 360 ymin ymax])
xlabel('angle'), ylabel('abs pattern')
