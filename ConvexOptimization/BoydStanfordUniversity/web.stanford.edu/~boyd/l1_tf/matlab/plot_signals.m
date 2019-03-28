function plot_signals(x,y,z1,z2,stockname,lambda,i,dates)
%
% plot_signals(x,y,z1,z2,stockname,lambda,i,dates)
%
% plot_signals displays original signal and fitted signals
%
%   INPUT
%       x           : n-vector; time index
%       y           : n-vector; original signal (time-series data)
%       z1          : n-vector; first fit
%       z2          : n-vector; second fit
%       stockname   : string; stock name to be used in figure file
%       lambda      : scalar; lambda vaule to be used in figure file
%       i           : scalar; figure start index (usually 1)
%       dates       : nx3 matrix; date matrix (year, month, day)
%

SAVE_FIGURE = false;

if (~SAVE_FIGURE)

    % trend plot
    figure(i);
    z = z1;
    subplot(3,1,1);
    plot(x,y,'k-','LineWidth',1.0); hold on;
    xlabel('date'); ylabel('log price');

    z = z1;
    subplot(3,1,2);
    plot(x,y,'k:','LineWidth',1.0); hold on;
    plot(x,z,'b-','LineWidth',2.0); hold off;
    xlabel('date'); ylabel('log price');

    z = z2;
    subplot(3,1,3);
    plot(x,y,'k:','LineWidth',1.0); hold on;
    plot(x,z,'b-','LineWidth',2.0); hold off;
    xlabel('date'); ylabel('log price');

else

    % date conversion
    idx = [false; diff(dates(:,1))~=0];
    xtick = x(idx);
    xticklabel = int2str(dates(idx,1));
    ytick = [floor(10*(min(y))):ceil(10*(max(y)))]'/10;
    yticklabel = num2str(ytick,'%2.1f');

    figure(i);
    z = z1;
    plot(x,y,'k-','LineWidth',1.0); hold on;
    xlabel('date'); ylabel('log price');
    set(gca,'FontSize',10,...
        'XTick', xtick, 'XTickLabel', xticklabel, ...
        'YTick', ytick, 'YTickLabel', yticklabel);
    print(i, '-depsc', sprintf('%s_org_%04d',stockname,floor(log10(lambda)*100)));


    figure(i+1);
    z = z1;
    plot(x,y,'k:','LineWidth',1.0); hold on;
    plot(x,z,'b-','LineWidth',2.0); hold off;
    xlabel('date'); ylabel('log price');

    set(gca,'FontSize',10,...
        'XTick', xtick, 'XTickLabel', xticklabel, ...
        'YTick', ytick, 'YTickLabel', yticklabel);
    print(i+1, '-depsc', sprintf('%s_l1tf_%04d',stockname,floor(log10(lambda)*100)));

    figure(i+2);
    z = z2;
    plot(x,y,'k:','LineWidth',1.0); hold on;
    plot(x,z,'b-','LineWidth',2.0); hold off;
    xlabel('date'); ylabel('log price');

    set(gca,'FontSize',10,...
        'XTick', xtick, 'XTickLabel', xticklabel, ...
        'YTick', ytick, 'YTickLabel', yticklabel);
    print(i+2, '-depsc', sprintf('%s_hp_%04d',stockname,floor(log10(lambda)*100)));

end