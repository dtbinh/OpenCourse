% Designs nominal and robust Chebychev FIR equalizer for
% a single-input single-output (SISO) channel as described
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
% can be formulated as a second-order cone program (SOCP).
%
% This script solves a numerical instance of the problem with
% a general ellipsoidal uncertainty (SDP problem).
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
% number of freq samples (rule-of-thumb)
m  = 10*(length(g) + n);

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
if strcmp(cvx_status,'Failed')
  return
end

% frequency response of the nominal equalizer
Hnom = A*hnom;

% nominal optimal value for the nominal equalizer
nom_nominal_optval = max( abs( Hnom.*G - Gdes ) );


%********************************************************************
% plots
%********************************************************************
% plot the channel freqency response (magnitude and phase)
figure(1), clf
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

% equalized frequency responses
Gnom_nom = Hnom.*G;
figure(3), clf
% magnitude
subplot(2,1,1);
plot(w,abs(Gnom_nom),'b-','LineWidth',1.5), hold on
plot(w,abs(Gdes),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
xlabel('w')
ylabel('mag H in dB')
axis([0 pi 0 2])
% phase
subplot(2,1,2)
plot(w,angle(Gnom_nom),'b-','LineWidth',1.5), hold on
plot(w,angle(Gdes),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
axis([0,pi,-pi,pi])
xlabel('w'), ylabel('phase H')

%********************************************************************
% pick a specific rho value and plot freq. responses
%********************************************************************
% amount of the multiplicative noise
rho = 0.05;
fprintf(1,'Solving specific robust equalization for rho = %3.4f\n',rho);

% solve robust equalization problem for the rho below
cvx_begin
  variable t
  variable hrob(n,1)

  minimize( t )
  subject to
    for k = 1:length(A);
      fprintf(1,'   * processing constraint %d\n',k);
      p = G(k); p = rho*norm(p)*p/norm(p);
      q = G(k); q = 4*rho*norm(q)*q/norm(q)*exp(i*pi/2);
      sup_ellip_equalizer( hrob,A(k,:),Gdes(k),G(k),p,q ) <= t;
    end
cvx_end

% frequency response of the robust equalizer
Hrob = A*hrob;

% plot nominal and robust equalizers
figure(2), clf
plot([0:n-1],hnom,'ro',[0:n-1],hnom,'r--'), hold on
plot([0:n-1],hrob,'bo',[0:n-1],hrob,'b-'), hold off
xlabel('n')
ylabel('h(n)')
set(gca,'FontSize',14);

% equalized frequency responses
Gnom_nom = Hnom.*G;
Grob_nom = Hrob.*G;
figure(3), clf
% magnitude
subplot(2,1,1);
plot(w,abs(Gnom_nom),'r--',w,abs(Grob_nom),'b-','LineWidth',1.5), hold on
plot(w,abs(Gdes),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
xlabel('w')
ylabel('mag H in dB')
axis([0 pi 0 2])
% phase
subplot(2,1,2)
plot(w,angle(Gnom_nom),'r--',w,angle(Grob_nom),'b-','LineWidth',1.5), hold on
plot(w,angle(Gdes),'-.','Color',[0 .5 0],'LineWidth',1.5), hold off
axis([0,pi,-pi,pi])
xlabel('w'), ylabel('phase H')
