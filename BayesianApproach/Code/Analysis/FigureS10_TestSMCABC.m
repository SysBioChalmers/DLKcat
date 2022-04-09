% This funciton is to generate the figure for Figure S4b
cd ../../Results/Test_SMCABC
nfound = length(dir('kcat_genra*.txt'));
for i = 1:nfound
    tmp = readmatrix(['kcat_genra',num2str(i),'.txt'],'FileType','text','Delimiter',',');
    theta_100(:,i) = tmp(end-1,:);
    theta_100_test(:,i) = tmp(end,:);
end

theta_meadin = median(theta_100);
tmp = sort(theta_100);
theta_5th = tmp(5,:);
theta_95th = tmp(95,:);

theta_test_median = median(theta_100_test);
tmp = sort(theta_100_test);
theta_test_5th = tmp(5,:);
theta_test_95th = tmp(95,:);

color = [33, 102, 172;178, 24, 43]./255;
hold on
x = [1:1:length(theta_test_median)]
plot(x,theta_test_median,'-','LineWidth',0.75,'Color',color(1,:))
X_plot = [x, fliplr(x)];
Y_plot = [theta_test_5th, fliplr(theta_test_95th)];
%nanx = isnan(X_plot)| isnan(Y_plot)
%X_plot = X_plot(~nanx)
%Y_plot = Y_plot(~nanx)

fill(X_plot, Y_plot , 1, ...
    'edgecolor','none', ...
    'facecolor',color(1,:), ...
    'edgecolor','none', ...
    'facealpha', 0.3);
plot(x,theta_meadin,'-','LineWidth',0.75,'Color',color(2,:))
X_plot = [x, fliplr(x)];
Y_plot = [theta_5th, fliplr(theta_95th)];
%nanx = isnan(X_plot)| isnan(Y_plot)
%X_plot = X_plot(~nanx)
%Y_plot = Y_plot(~nanx)
Y_plot(isnan(Y_plot)) = 0;
fill(X_plot, Y_plot , 1, ...
    'edgecolor','none', ...
    'facecolor',color(2,:), ...
    'edgecolor','none', ...
    'facealpha', 0.3);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('RMSE in each generation','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('No. generation','FontSize',7,'FontName','Helvetica','Color','k');
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);
box off
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
'XAxisLocation','top', 'XTick', [],...
'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])

cd ../../Code/Analysis