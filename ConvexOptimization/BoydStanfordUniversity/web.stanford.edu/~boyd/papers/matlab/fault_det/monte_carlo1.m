% MONTE_CARLO1.M
%
% Generates Figure 3 from
% "Relaxed Maximum a Posteriori Fault Identification"
% by A. Zymnis, S. Boyd and D. Gorinevsky
%
% Note: Number of samples is large, might take a while to run

clear all

randn('state',0);
rand('state',0);

%% ---------------- Generate Data -----------------------------------------
m = 50; %number of sensors
n = 100; %number of fault signatures
K = 10;
pf = 0.05; %probability of fault
sigmas = linspace(0.1,4,30); %noise std
lambda = log((1-pf)/pf)*ones(n,1);

%% --------------- Generate Monte Carlo Examples ---------------------
N_samples = 2000;
S1_sort = []; S10_sort = [];
S1_lc = []; S10_lc = [];
iter_sim = 1;
for j = 1:length(sigmas)
    sigma = sigmas(j);
    success1_sort = 0; success10_sort = 0;
    success1_lc = 0; success10_lc = 0;
    for iter = 1:N_samples

        A = randn(m,n); %fault signatures
        x_true = double(rand(n,1)<pf); %true fault vector
        X_true = x_true*ones(1,K);
        y = A*x_true+sigma*randn(m,1); %measurement
        
        %% Without local optimization
        [X_amb,l_amb] = armap(A,y,pf*ones(n,1),sigma,1/(2*n),0,'quiet');
        s1a = all(x_true==X_amb(:,1));
        success1_sort = s1a+success1_sort;
        s10a = any(all(X_true==X_amb));
        success10_sort = s10a+success10_sort;
        
        %% With local optimization
        [X_amb_lc,l_amb_lc] = armap(A,y,pf*ones(n,1),sigma,1/(2*n),1,'quiet');
        s1lc = all(x_true==X_amb_lc(:,1));
        s10lc = any(all(X_true==X_amb_lc));
        success1_lc = success1_lc+s1lc;
        success10_lc = success10_lc+s10lc;

        fprintf(1,'Iter %3d of %d...sigma = %3.3f, SORT: %d %d, LC: %d %d\n',...
            iter_sim,N_samples*length(sigmas),sigma,s1a,s10a,s1lc,s10lc)
        iter_sim = iter_sim+1;
    end
    S1_sort = [S1_sort 100*success1_sort/N_samples];
    S10_sort = [S10_sort 100*success10_sort/N_samples];
    S1_lc = [S1_lc 100*success1_lc/N_samples];
    S10_lc = [S10_lc 100*success10_lc/N_samples];
    save sigma_comp_lc sigmas S1_sort S10_sort S1_lc S10_lc
end
save sigma_comp_lc sigmas S1_sort S10_sort S1_lc S10_lc

figure
subplot(2,1,1)
semilogx(sqrt(pf*n)./sigmas,S1_sort,'-',sqrt(pf*n)./sigmas,S10_sort,'--','linewidth',2)
axis([0.6 5 0 100])
set(gca,'xtick',[0.6 1 2 3 4 5])
set(gca,'fontsize',12)
title('SRMAP')
ylabel('probsuccess')
subplot(2,1,2)
semilogx(sqrt(pf*n)./sigmas,S1_lc,'-',sqrt(pf*n)./sigmas,S10_lc,'--','linewidth',2)
axis([0.6 5 0 100])
set(gca,'xtick',[0.6 1 2 3 4 5])
set(gca,'fontsize',12)
title('locopt')
xlabel('SNR')
ylabel('probsuccess')
%print -depsc sigma_res.eps

figure
semilogx(sqrt(pf*n)./sigmas,S1_sort,'-',sqrt(pf*n)./sigmas,S1_lc,'--','linewidth',2)
axis([0.6 5 0 100])
set(gca,'xtick',[0.6 1 2 3 4 5])
set(gca,'fontsize',12)
xlabel('SNR')
ylabel('probsuccess')
%print -depsc sigma_1_lc_comp.eps

figure
semilogx(sqrt(pf*n)./sigmas,S10_sort,'-',sqrt(pf*n)./sigmas,S10_lc,'--','linewidth',2)
axis([0.6 5 0 100])
set(gca,'xtick',[0.6 1 2 3 4 5])
set(gca,'fontsize',12)
xlabel('SNR')
ylabel('probsuccess')
%print -depsc sigma_10_lc_comp.eps