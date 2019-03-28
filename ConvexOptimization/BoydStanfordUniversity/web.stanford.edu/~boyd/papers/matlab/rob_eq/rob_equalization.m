% The following Matlab script generates all the figures 
% in the paper: "Robust Chebychev FIR equalization"
% by A. Mutapcic, S.-J. Kim, and S. Boyd, Proc. Globecom 07.
%
% The nominal equalizer is a solution of the convex optimization problem:
%
%   minimize   max_m |G(w_m)H(w_m) - G_des(w_m)|  for w_m in [0,pi] 
%
% where the optimization variable is the frequency response
% function H and therefore the filter impulse response h.
% The problem data is the nominal unequalized frequency
% response G and the desired frequency response G_des.
% The frequency is sampled at M points w_1, ... , w_M.
%
% In robust equalization, we consider the case when G(w_m) is uncertain,
% but known to lie in an ellipsoidal uncertainty set E_m for each w_m.
% We solve the following (worst-case) robust equalization problem:
%
%   minimize   max_m sup_{E_m} |G(w_m)H(w_m) - G_des(w_m)|  for w_m in [0,pi] 
%
% The above problem can be cast as an SDP (see the paper).
% When the uncertainty is given by a complex disk, the design problem
% can be formulated as a second-order cone program (SOCP);
% this script solves a numerical instance of such a problem (SOCP).
%
% Almir Mutapcic, Sept 2007.

%********************************************************************
% problem specs
%********************************************************************
% sample channel with impulse response g
g = 0.5*[0.25; .5; 0; 1; 0.25;];

% problem parameters
n   = 20;   % filter order
D   = 8;    % overall delay

%********************************************************************
% problem data structures
%********************************************************************
% number of freq samples
m = 100;

w = linspace(0,pi,m)';
A = exp( -i*kron(w,[0:n-1]) );           % DFT matrix
G = exp( -i*kron(w,[0:length(g)-1]) )*g; % nominal freq response

% desired frequency response is a pure delay (equalized channel)
Gdes = exp(-i*D*w);

%********************************************************************
% nominal equalization (ignores uncertainty)
%********************************************************************
% formulate and solve the Chebyshev design problem
cvx_begin
  variable hnom(n,1)
  minimize( max( abs( (A*hnom).*G - Gdes ) ) ) 
cvx_end

% check if problem was successfully solved
disp(['Frequency equalization problem is ' cvx_status])
if ~strcmp(cvx_status,'Solved')
  return
end

% frequency response of the nominal equalizer
Hnom = A*hnom;

% nominal optimal value for the nominal equalizer
nom_nominal_optval = max( abs( Hnom.*G - Gdes ) );

%********************************************************************
% robust equalization (takes into account uncertainty)
%********************************************************************
% solve the robust problem for rho varying between 0 and 0.1
rho = linspace(0,0.1,20);
invG = 1./abs(G);
rho_pattern = invG/max(invG);

nom_robust_optval = [];
wc_robust_optval  = [];
wc_nominal_optval = [];

cvx_quiet(true);
for ind = 1:length(rho)
  fprintf(1,'  *** Solving robust equalization for rho = %3.4f\n',rho(ind));
  rho_vect = rho(ind)*rho_pattern;
  cvx_begin
    variable hrob(n,1)
    minimize( max( abs( (A*hrob).*G - Gdes ) + rho_vect.*abs( A*hrob ) ) ) 
  cvx_end

  % frequency response of the robust equalizer
  Hrob = A*hrob;

  % collect performance data for the nominal and the robust equalizers
  nom_robust_optval(end+1) = max( abs( Hrob.*G - Gdes ) );
  wc_robust_optval(end+1)  = max( abs( Hrob.*G - Gdes ) + rho_vect.*abs( Hrob ) );
  wc_nominal_optval(end+1) = max( abs( Hnom.*G - Gdes ) + rho_vect.*abs( Hnom ) );
end
cvx_quiet(false);

%********************************************************************
% plots
%********************************************************************
% plot the channel impulse response
figure(1), clf
stem([0:length(g)-1],g)
axis([-.5 4.5 0 0.75])
xlabel('n')
ylabel('g(n)')

% plot the channel freqency response (magnitude and phase)
figure(2), clf
% magnitude
subplot(2,1,1);
plot(w,20*log10(abs(G)),'k-','LineWidth',1.5), hold on
plot(w,20*log10(abs(Gdes)),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
axis([0,pi,-15,5])
ylabel('mag G(w) in dB')
% phase
subplot(2,1,2)
plot(w,angle(G),'k-','LineWidth',1.5), hold on
plot(w,angle(Gdes),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
axis([0,pi,-pi,pi])
xlabel('w')
ylabel('phase G(w)')

% worst-case rejection level versus rho
figure(3), clf
plot(rho,wc_robust_optval,'LineWidth',1.5), hold on
plot(rho,wc_nominal_optval,'r--','LineWidth',1.5), hold off,
xlabel('rho'), ylabel('worst-case MAD')
axis([0 0.1 0.1 0.6]), grid on

figure(4), clf
plot(rho,nom_robust_optval,'LineWidth',1.5), hold on
plot(rho,nom_nominal_optval*ones(length(rho),1),'r--','LineWidth',1.5), hold off,
xlabel('rho'), ylabel('nominal objective value')
axis([0 0.1 0.1 0.6]), grid on

%********************************************************************
% pick a specific rho value and plot freq. responses
%********************************************************************
rho = 0.1;
rho_vect = rho*rho_pattern;
fprintf(1,'Solving specific robust equalization for rho = %3.4f\n',rho);
% solve robust equalization problem for the rho below
cvx_begin
  variable hrob(n,1)
  minimize( max( abs( (A*hrob).*G - Gdes ) + rho_vect.*abs( A*hrob ) ) ) 
cvx_end

% frequency response of the robust equalizer
Hrob = A*hrob;

% plot nominal and robust equalizers
figure(5), clf
plot([0:n-1],hnom,'ro',[0:n-1],hnom,'r--'), hold on
plot([0:n-1],hrob,'bo',[0:n-1],hrob,'b-'), hold off
xlabel('n')
ylabel('h(n)')

% equalized frequency responses
Gnom_nom = Hnom.*G;
Grob_nom = Hrob.*G;
figure(6), clf
% magnitude
subplot(2,1,1);
plot(w,20*log10(abs(Gnom_nom)),'r--',w,20*log10(abs(Grob_nom)),'b-','LineWidth',1.5), hold on
plot(w,20*log10(abs(Gdes)),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
xlabel('w')
ylabel('mag H in dB')
axis([0 pi -15 5])
% phase
subplot(2,1,2)
plot(w,angle(Gnom_nom),'r--',w,angle(Grob_nom),'b-','LineWidth',1.5), hold on
plot(w,angle(Gdes),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
axis([0,pi,-pi,pi])
xlabel('w'), ylabel('phase H')

% plot the worst-case frequency responses
Gnom_wc = ( G + rho_vect.*exp(i*angle(Hnom.*G-Gdes)-i*angle(Hnom)) ).*Hnom;
Grob_wc = ( G + rho_vect.*exp(i*angle(Hrob.*G-Gdes)-i*angle(Hrob)) ).*Hrob;
figure(7), clf
% magnitude
subplot(2,1,1);
plot(w,20*log10(abs(Gnom_wc)),'r--',w,20*log10(abs(Grob_wc)),'b-','LineWidth',1.5), hold on
plot(w,20*log10(abs(Gdes)),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
xlabel('w')
ylabel('mag H in dB')
axis([0 pi -15 5])
% phase
subplot(2,1,2)
plot(w,angle(Gnom_wc),'r--',w,angle(Grob_wc),'b-','LineWidth',1.5), hold on
plot(w,angle(Gdes),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
axis([0,pi,-pi,pi])
xlabel('w'), ylabel('phase H')
