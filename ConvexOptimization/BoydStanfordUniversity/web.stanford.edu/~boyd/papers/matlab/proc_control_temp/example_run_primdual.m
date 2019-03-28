% EXAMPLE_RUN_PRIDUAL  solves the processor speed problem using the primal-dual method from
%                      "Processor speed control with thermal constraints"
%                      by Mutapcic, Boyd, Murali, Atienza, De Micheli, Gupta
%
% Requires CVX for optimal solution comparison.

% almirm, 9/07

% create problem data for the example in the paper
create_example;

%
% primal-dual solution
%
% call our custom primal-dual algorithm
fprintf(1,'Calling custom primal-dual interior-point method ...\n');
[spd,lamspd,iters_pd,gaps_pd,residl_pd] = pd_thr_speed(G,phi,Tother,Tamb,Tmax,smin,smax,1e-8,0);

%
% use CVX to find solution and check the optimal values
%
fprintf(1,'\nCalling CVX to solve the same problem, so we can check answers ...\n');
cvx_begin
  variable s_opt(n)    % processor speed variables
  maximize sum(s_opt)  % maximize throughput
  subject to
    % maximum temperature constraint
    G*pow_pos(s_opt,3) + Tother + Tamb <= Tmax;
    % processor speed bounds
    smin <= s_opt; s_opt <= smax;
cvx_end
optval = cvx_optval;
fprintf(1,'Optimal maximum workload is %3.4f\n',optval);

fprintf(1,'\n\n* abs  error between optimal values    is %3.4f\n',abs(sum(spd)-optval));
fprintf(1,'* norm error between optimal solutions is %3.4f\n\n',norm(spd-s_opt));

%
% suboptimal solution (find optimal constant speed for all processors)
%
[sub_optval, s_const] = opt_equal_speed(G,Tother,Tamb,Tmax,smin,smax);;
fprintf(1,'Sub-optimal maximum workload is %3.4f with constant s = %3.4f\n',...
           sub_optval,s_const);

%
% plot results 
%
figure(1)
iter = [1:iters_pd];
semilogy(iter,gaps_pd,'r-o',iter,residl_pd,'b--s','LineWidth',1.15);
xlabel('iters');
legend('duality measure','residualn')

% optimal solution
figure(2), clf
temps = reshape(G*s_opt.^3 + Tother + Tamb, node_m, node_n);
imagesc(temps,[25 75])
colorbar
axis([1 node_n 1 node_m]),
set(gca,'XTick',[5.5:5:node_n],'YTick',[5.5:5:node_m])
set(gca,'XTickLabel',[],'YTickLabel',[])
title(sprintf('optval %3.2f',optval))

% suboptimal solution
figure(3), clf
temps = reshape(G*ones(n,1)*s_const^3 + Tother + Tamb, node_m, node_n);
imagesc(temps,[25 75])
colorbar
axis([1 node_n 1 node_m]),
set(gca,'XTick',[5.5:5:node_n],'YTick',[5.5:5:node_m])
set(gca,'XTickLabel',[],'YTickLabel',[])
title(sprintf('optval %3.2f',sub_optval))

% histogram of optimal solution frequencies (vs. optimal constant freq)
figure(4), clf
edges = linspace(0.75,3.25,15);
[count] = histc( s_opt, edges )';
bar( edges, count, 'histc' );
hold on
h = findobj(gca,'Type','line'); delete(h)
hpatch = findobj(gca,'Type','patch');
set(hpatch,'FaceColor',[0.85,0.85,0.85])
set(gca,'FontSize',13);
axis([0.75 3.25 0 17.5])
set(gca,'YTick',[0 5 10 15]);
plot([s_const s_const],[0 17.5],'k','LineWidth',2),
text(s_const+.05,14,'fconst')
hold off
