% EXAMPLE_RUN_DISTR  solves the processor speed problem using distributed method from
%                    "Processor speed control with thermal constraints"
%                    by Mutapcic, Boyd, Murali, Atienza, De Micheli, Gupta

% almirm, 9/07

% create problem data for the example in the paper
create_example;

% set seed number so we can reproduce results
rand('state',0);

% two scenarios:
% the first one is disturbed version of nominal Tother (corresponds
%               to about 10% distrubance in both unknown power sources
%               and about 10% rise in Tamb -- this is folded into Tother)
% the second one is the nominal Tother from the create_example script
Tother1 = Tother + 5*rand(m,1);
Tother2 = Tother;

%
% compute the optimal values for both scenarios using the primal-dual algorithm
%
fprintf(1,'Calling custom primal-dual interior-point method for scenario 1...\n');
[spd,lamspd,iters_pd,gaps_pd,residl_pd] = pd_thr_speed(G,phi,Tother1,Tamb,Tmax,smin,smax,1e-8,0);
optval1 = sum(spd);

fprintf(1,'\nCalling custom primal-dual interior-point method for scenario 2...\n');
[spd,lamspd,iters_pd,gaps_pd,residl_pd] = pd_thr_speed(G,phi,Tother2,Tamb,Tmax,smin,smax,1e-8,0);
optval2 = sum(spd);;


%
% simple gradient method for the dual problem
%
NUM_ITERS1 = 200;         % scenario switch occurs at k = 200
NUM_ITERS2 = 400;
NUM_ITERS  = NUM_ITERS2;

alpha  = 1e-4;
lambda = 0.1*ones(m,1);

primal_hist  = zeros(NUM_ITERS,1);
maxtemp_hist = zeros(NUM_ITERS,1);
norm_hist    = zeros(NUM_ITERS,1);

% first scenario
for iter = 1:NUM_ITERS1
  % individual updates
  fprintf(1,'Iteration number %d\n',iter)
  s_dist = 1./sqrt(3*G'*lambda);
  s_dist = max( min(s_dist, smax), smin);
  T = G*s_dist.^3 + Tother1 + Tamb;
  lambda = max( lambda - alpha*(Tmax-T), 0 );

  % book-keeping
  primal_hist(iter)  = sum(s_dist);
  maxtemp_hist(iter) = max(T);
  norm_hist(iter)    = norm( s_opt - s_dist );
end

% second scenario
for iter = NUM_ITERS1+1:NUM_ITERS2
  % individual updates
  fprintf(1,'Iteration number %d\n',iter)
  s_dist = 1./sqrt(3*G'*lambda);
  s_dist = max( min(s_dist, smax), smin);
  T = G*s_dist.^3 + Tother2 + Tamb;
  lambda = max( lambda - alpha*(Tmax-T), 0 );

  % book-keeping
  primal_hist(iter)  = sum(s_dist);
  maxtemp_hist(iter) = max(T);
  norm_hist(iter)    = norm( s_opt - s_dist );
end

iters = [1:NUM_ITERS2];
% plot primal objective and maximum temperature on the chip at each iteration
figure(1), clf
subplot(2,1,1),
plot(iters, primal_hist(1:length(iters)),'LineWidth',1.5), hold on
plot([1:NUM_ITERS1], optval1*ones(NUM_ITERS1,1),'--','Color',[0 0.5 0],'LineWidth',1.5),
plot([NUM_ITERS1+1:NUM_ITERS2], optval2*ones(NUM_ITERS2-NUM_ITERS1,1),...
      '--','Color',[0 0.5 0],'LineWidth',1.25),
hold off
text(325,optval2-5,'Uopt')
ylabel('Uf')
axis([1 NUM_ITERS 150 190])
set(gca,'XTick',[0:100:400])

subplot(2,1,2)
plot(iters, maxtemp_hist(1:length(iters)),'LineWidth',1.25), hold on
plot(iters, Tmax*ones(length(iters),1),'r--','LineWidth',1.25), hold off
text(325,Tmax+3,'Tmax')
xlabel('iters'), ylabel('max temp')
axis([1 NUM_ITERS 60 90])
set(gca,'XTick',[0:100:400])


%
% smoothing gradient method for the dual problem
%

alpha  = 0.03;
theta  = alpha;
lambda = ones(m,1);
s_prev = smin*ones(n,1);

% first scenario
for iter = 1:NUM_ITERS1
  % individual updates
  fprintf(1,'Iteration number %d\n',iter)
  s_dist = 1./sqrt(3*G'*lambda);
  s_dist = max( min(s_dist, smax), smin);
  s_dist = theta*s_dist + (1-theta)*s_prev;
  s_prev = s_dist;
  T = G*s_dist.^3 + Tother1 + Tamb;
  lambda = max( lambda - alpha*(Tmax-T), 0 );

  % book-keeping
  primal_hist(iter)  = sum(s_dist);
  maxtemp_hist(iter) = max(T);
  norm_hist(iter)    = norm( s_opt - s_dist );
end

% second scenario
for iter = NUM_ITERS1+1:NUM_ITERS2
  % individual updates
  fprintf(1,'Iteration number %d\n',iter)
  s_dist = 1./sqrt(3*G'*lambda);
  s_dist = max( min(s_dist, smax), smin);
  s_dist = theta*s_dist + (1-theta)*s_prev;
  s_prev = s_dist;
  T = G*s_dist.^3 + Tother2 + Tamb;
  lambda = max( lambda - alpha*(Tmax-T), 0 );

  % book-keeping
  primal_hist(iter)  = sum(s_dist);
  maxtemp_hist(iter) = max(T);
  norm_hist(iter)    = norm( s_opt - s_dist );
end

iters = [1:NUM_ITERS2];
% plot primal objective and maximum temperature on the chip at each iteration
figure(2), clf
subplot(2,1,1),
plot(iters, primal_hist(1:length(iters)),'LineWidth',1.5), hold on
plot([1:NUM_ITERS1], optval1*ones(NUM_ITERS1,1),'--','Color',[0 0.5 0],'LineWidth',1.5),
plot([NUM_ITERS1+1:NUM_ITERS2], optval2*ones(NUM_ITERS2-NUM_ITERS1,1),...
      '--','Color',[0 0.5 0],'LineWidth',1.25),
hold off
text(325,optval2-5,'Uopt')
ylabel('Uf')
axis([1 NUM_ITERS 150 190])
set(gca,'XTick',[0:100:400])

subplot(2,1,2)
plot(iters, maxtemp_hist(1:length(iters)),'LineWidth',1.25), hold on
plot(iters, Tmax*ones(length(iters),1),'r--','LineWidth',1.25), hold off
text(325,Tmax+3,'Tmax')
xlabel('iters'), ylabel('max temp')
axis([1 NUM_ITERS 60 90])
set(gca,'XTick',[0:100:400])
