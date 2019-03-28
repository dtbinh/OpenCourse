% MONTE_CARLO_SDP.M
%
% Generates Figure 3 from
% "Relaxed Maximum a Posteriori Fault Identification"
% by A. Zymnis, S. Boyd and D. Gorinevsky
%
% Note: Number of samples is large, might take a while to run

clear all

randn('state',0);
rand('state',0);
cvx_quiet(true);

%% ---------------- Generate Data -----------------------------------------
m = 50; %number of sensors
n = 100; %number of fault signatures
K = 10;
pf = 0.05; %probability of fault
sigmas = linspace(0.1,4,30); %noise std
lambda = log((1-pf)/pf)*ones(n,1);

%% --------------- Generate Monte Carlo Examples ---------------------
N_samples = 2000;
S1_sdp = []; S10_sdp = [];
S1_sdp_lc = []; S10_sdp_lc = [];
iter_sim = 1;
for j = 1:length(sigmas)
    sigma = sigmas(j);
    success1_sdp = 0; success10_sdp = 0;
    success1_sdp_lc = 0; success10_sdp_lc = 0;
    for iter = 1:N_samples

        A = randn(m,n); %fault signatures
        x_true = double(rand(n,1)<pf); %true fault vector
        X_true = x_true*ones(1,K);
        y = A*x_true+sigma*randn(m,1); %measurement
        
        % solve SDP relaxation
        W = (1/(2*sigma^2))*A'*A;
        w = (lambda-(1/sigma^2)*A'*y);
        cvx_begin
            variable x(n)
            variable X(n,n) symmetric
            minimize(trace(W*X)+w'*x)
            subject to
                diag(X)==x;
                [X x; x' 1] == semidefinite(n+1);
                x >= 0;
                x <= 1;
        cvx_end
        l_min_sdp = cvx_optval;

        % rounding
        [x_sort,ind_x] = sort(x,'descend'); x_cand = []; l_cand = [];
        for i = 1:n
            x_cur = zeros(n,1);
            x_cur(ind_x(1:i)) = 1;
            l_cur = (1/(2*sigma^2))*square_pos(norm(A*x_cur-y,2))...
            +lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
            x_cand = [x_cand x_cur]; l_cand = [l_cand l_cur]; 
        end
        [l_sort,ind_l] = sort(l_cand,'ascend');
        X_amb_sdp = x_cand(:,ind_l(1:K)); %get ambiguity set
        l_amb_sdp = l_sort(1:K);


        s1a = all(x_true==X_amb_sdp(:,1));
        success1_sdp = s1a+success1_sdp;
        s10a = any(all(X_true==X_amb_sdp));
        success10_sdp = s10a+success10_sdp;


        % perform local optimization
        EXIT_FLAG = 0; iter = 0;
        while(~EXIT_FLAG)
            x_cur = X_amb_sdp(:,1); x_best = x_cur;
            for i = 1:n
                iter = iter+1;
                x_cur(i) = not(x_cur(i));
                l_cur = (1/(2*sigma^2))*norm(A*x_cur-y,2).^2+lambda'*x_cur-(1/(2*sigma^2))*norm(y)^2;
                if any(l_cur<l_amb_sdp)
                    ind = find(l_cur<l_amb_sdp);
                    ind = ind(1);
                    X_amb_sdp = [X_amb_sdp(:,1:(ind-1)) x_cur X_amb_sdp(:,ind:(end-1))];
                    l_amb_sdp = [l_amb_sdp(:,1:(ind-1)) l_cur l_amb_sdp(:,ind:(end-1))];
                    if ind==1
                        fprintf(1,'Found new best pattern!\n');
                    else
                        x_cur(i) = not(x_cur(i));
                    end
                else
                    x_cur(i) = not(x_cur(i));
                end
            end
            if all(x_best == X_amb_sdp(:,1)), 
                EXIT_FLAG = 1; 
            end
        end

        s1lc = all(x_true==X_amb_sdp(:,1));
        success1_sdp_lc = s1lc+success1_sdp_lc;
        s10lc = any(all(X_true==X_amb_sdp));
        success10_sdp_lc = s10lc+success10_sdp_lc;
 
        fprintf(1,'Iter %3d of %d, SDP...sigma = %3.3f, SORT: %d %d, LC: %d %d\n',...
            iter_sim,N_samples*length(sigmas),sigma,s1a,s10a,s1lc,s10lc)
        iter_sim = iter_sim+1;
    end
    S1_sdp = [S1_sdp 100*success1_sdp/N_samples];
    S10_sdp = [S10_sdp 100*success10_sdp/N_samples];
    S1_sdp_lc = [S1_sdp_lc 100*success1_sdp_lc/N_samples];
    S10_sdp_lc = [S10_sdp_lc 100*success10_sdp_lc/N_samples];
    save sigma_sdp sigmas S1_sdp S10_sdp S1_sdp_lc S10_sdp_lc
end
save sigma_sdp sigmas S1_sdp S10_sdp S1_sdp_lc S10_sdp_lc

return

figure
subplot(2,1,1)
semilogx(sqrt(pf*n)./sigmas,S1_sdp,'-',sqrt(pf*n)./sigmas,S10_sdp,'--','linewidth',2)
axis([0.6 5 0 100])
set(gca,'xtick',[0.6 1 2 3 4 5])
set(gca,'fontsize',12)
title('SRMAP')
ylabel('probsuccess')
subplot(2,1,2)
semilogx(sqrt(pf*n)./sigmas,S1_sdp_lc,'-',sqrt(pf*n)./sigmas,S10_sdp_lc,'--','linewidth',2)
axis([0.6 5 0 100])
set(gca,'xtick',[0.6 1 2 3 4 5])
set(gca,'fontsize',12)
title('locopt')
xlabel('SNR')
ylabel('probsuccess')
%print -depsc sigma_res.eps

figure
semilogx(sqrt(pf*n)./sigmas,S1_sdp,'-',sqrt(pf*n)./sigmas,S1_sdp_lc,'--','linewidth',2)
axis([0.6 5 0 100])
set(gca,'xtick',[0.6 1 2 3 4 5])
set(gca,'fontsize',12)
xlabel('SNR')
ylabel('probsuccess')
%print -depsc sigma_1_lc_comp.eps

figure
semilogx(sqrt(pf*n)./sigmas,S10_sdp,'-',sqrt(pf*n)./sigmas,S10_sdp_lc,'--','linewidth',2)
axis([0.6 5 0 100])
set(gca,'xtick',[0.6 1 2 3 4 5])
set(gca,'fontsize',12)
xlabel('SNR')
ylabel('probsuccess')
%print -depsc sigma_10_lc_comp.eps
