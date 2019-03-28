% MONTE_CARLO2.M
%
% Generates Figure 4 from
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

%% --------------- Generate 1000 Monte Carlo Examples ---------------------
N_samples = 4000;
S1_double = []; S10_double = [];
S1_1bit = []; S10_1bit = [];
S1_2bit = []; S10_2bit = [];
S1_3bit = []; S10_3bit = [];
S1_4bit = []; S10_4bit = [];
iter_sim = 1;
for j = 1:length(sigmas)
    success1_double = 0; success10_double = 0;
    success1_1bit = 0; success10_1bit = 0;
    success1_2bit = 0; success10_2bit = 0;
    success1_3bit = 0; success10_3bit = 0;
    success1_4bit = 0; success10_4bit = 0;
    completed1 = 0; completed2 = 0; completed3 = 0; completed4 = 0;
    for iter = 1:N_samples
        sigma = sigmas(j);
        A = randn(m,n); %fault signatures
        x_true = double(rand(n,1)<pf); %true fault vector
        X_true = x_true*ones(1,K);
        y = A*x_true+sigma*randn(m,1); %measurement
        
        
        %% DOUBLE
        [X_amb,l_amb,ERROR_FLAG] = armap(A,y,pf*ones(n,1),sigma,0.05,1,'quiet');
        
        s1d = all(x_true==X_amb(:,1));
        success1_double = s1d+success1_double;
        s10d = any(all(X_true==X_amb));
        success10_double = s10d+success10_double;
        
        %% 1BIT
        % Quantize measurement
        t = 0;
        interval = sum(double(y*ones(1,length(t))>ones(m,1)*t),2);
        t = [-10 t 10];
        l = t(interval+1)'; u = t(interval+2)'; %measurement intervals
        y1 = (min(u,5)+max(l,-5))/2;;
        
        [X_amb,l_amb,ERROR_FLAG] = armap(A,[l u],pf*ones(n,1),sigma,0.05,1,'quiet');
        
        s1_1 = 0; s10_1 = 0;
        if ~ERROR_FLAG
            completed1 = completed1+1;
            s1_1 = all(x_true==X_amb(:,1));
            success1_1bit = s1_1+success1_1bit;
            s10_1 = any(all(X_true==X_amb));
            success10_1bit = s10_1+success10_1bit;
        end
        
        %% 2BIT
        % Quantize measurement
        t = [-2 0 2];
        interval = sum(double(y*ones(1,length(t))>ones(m,1)*t),2);
        t = [-10 t 10];
        l = t(interval+1)'; u = t(interval+2)'; %measurement intervals
        y2 = (min(u,4)+max(l,-4))/2;;
        
        [X_amb,l_amb,ERROR_FLAG] = armap(A,[l u],pf*ones(n,1),sigma,0.05,1,'quiet');
        
        s1_2 = 0; s10_2 = 0;
        if ~ERROR_FLAG
            completed2 = completed2+1;
            s1_2 = all(x_true==X_amb(:,1));
            success1_2bit = s1_2+success1_2bit;
            s10_2 = any(all(X_true==X_amb));
            success10_2bit = s10_2+success10_2bit;
        end
        
        %% 3BIT
        % Quantize measurement
        t = [-6 -4 -2 0 2 4 6];
        interval = sum(double(y*ones(1,length(t))>ones(m,1)*t),2);
        t = [-10 t 10];
        l = t(interval+1)'; u = t(interval+2)'; %measurement intervals
        y3 = (min(u,8)+max(l,-8))/2;;
        
        [X_amb,l_amb,ERROR_FLAG] = armap(A,[l u],pf*ones(n,1),sigma,0.05,1,'quiet');
        
        s1_3 = 0; s10_3 = 0;
        if ~ERROR_FLAG
            completed3 = completed3+1;
            s1_3 = all(x_true==X_amb(:,1));
            success1_3bit = s1_3+success1_3bit;
            s10_3 = any(all(X_true==X_amb));
            success10_3bit = s10_3+success10_3bit;
        end
 
        %% 4BIT
        % Quantize measurement
        t = [-7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7];
        interval = sum(double(y*ones(1,length(t))>ones(m,1)*t),2);
        t = [-10 t 10];
        l = t(interval+1)'; u = t(interval+2)'; %measurement intervals
        y4 = (min(u,8.5)+max(l,-8.5))/2;
        
        [X_amb,l_amb,ERROR_FLAG] = armap(A,[l u],pf*ones(n,1),sigma,0.05,1,'quiet');
        
        s1_4 = 0; s10_4 = 0;
        if ~ERROR_FLAG
            completed4 = completed4+1;
            s1_4 = all(x_true==X_amb(:,1));
            success1_4bit = s1_4+success1_4bit;
            s10_4 = any(all(X_true==X_amb));
            success10_4bit = s10_4+success10_4bit;
        end
        
        fprintf(1,'Iter %3d of %d...sigma = %3.3f, DOUBLE: %d %d, 1BIT: %d %d, 2BIT: %d %d, 3BIT: %d %d, 4BIT: %d %d\n',...
            iter_sim,N_samples*length(sigmas),sigmas(j),s1d,s10d,s1_1,s10_1,s1_2,s10_2,s1_3,s10_3,s1_4,s10_4)
        iter_sim = iter_sim+1;
    end
    S1_double = [S1_double 100*success1_double/N_samples];
    S10_double = [S10_double 100*success10_double/N_samples];
    S1_1bit = [S1_1bit 100*success1_1bit/completed1];
    S10_1bit = [S10_1bit 100*success10_1bit/completed1];
    S1_2bit = [S1_2bit 100*success1_2bit/completed2];
    S10_2bit = [S10_2bit 100*success10_2bit/completed2];
    S1_3bit = [S1_3bit 100*success1_3bit/completed3];
    S10_3bit = [S10_3bit 100*success10_3bit/completed3];
    S1_4bit = [S1_4bit 100*success1_4bit/completed4];
    S10_4bit = [S10_4bit 100*success10_4bit/completed4];
end
save sigma_results_bits_lc sigmas S1_double S10_double ...
    S1_1bit S10_1bit S1_2bit S10_2bit S1_3bit S10_3bit S1_4bit S10_4bit 

return

figure
SNR = sqrt(pf*n)./sigmas;
semilogx(SNR,S1_double,SNR,S1_1bit,SNR,S1_2bit,SNR,S1_3bit,SNR,S1_4bit,'linewidth',2)
axis([0.6 5 0 100])
xlabel('SNR'); ylabel('probsuccess')
set(gca,'fontsize',15)
set(gca,'xtick',[0.6 1 2 3 4 5])
print -depsc bits_compare.eps


