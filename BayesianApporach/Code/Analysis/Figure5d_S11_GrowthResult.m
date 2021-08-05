% Figure 5d
% This function is to plot the figure5d, which is mainly focusing on
% evalulation of the kcat tuning process
% Running this function requires that all tuned emodels are genrated

% [~,~,growthrates] = xlsread('growthratedata.xlsx','growthrates');
% growthrates = growthrates(2:end,:);
% species = unique(growthrates(:,4));
% species = strrep(species,' ','_');
load('species_withdata.mat')
species = species_withdata;

% Figure 3A this plot is to show the growth rate prediction vs experimental
% values
result_constrain = [];
result_max = [];
cd ../../Results/model_build_files/model_Bayesian
for i = 1:length(species)
    cd(species{i})% later change to cd(species{i}) use species id to name the folder
    load([species{i},'_sim_phen.mat']) % is generated in the function simulategrowth.m
    result_constrain = [result_constrain;growthdata(:,[1:3,14]),num2cell(simulated(1:length(growthdata(:,1)),1))];
    result_max = [result_max;max_growth(:,[1:3,14]),num2cell(simulated(length(growthdata(:,1))+1:end,1))];
    cd ../
end

% datapoint
all_data = [result_constrain;result_max];
dp = cell2mat(all_data(:,[3,5]));

% Figure 5d
hold on
plot(dp(:,1),dp(:,2),'o','LineWidth',0.75,'MarkerSize',3,'Color',[197,27,138]./255);

set(gca,'FontSize',6,'FontName','Helvetica');
limit = max(dp(:));
xlim([0 limit]);
ylim([0 limit]);
x = [0:0.1:limit];
plot(x,x,'--k','LineWidth',0.75)
set(gca,'FontSize',6,'FontName','Helvetica');

ylabel('Simulated growth rate [1/h]','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('Experimental growth rate [1/h]','FontSize',7,'FontName','Helvetica','Color','k');

set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
'XAxisLocation','top', 'XTick', [],...
'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])


% Supplementary Figure S11 split by condition 
sub = unique(all_data(:,2));
color = [166,206,227
    31,120,180
178,223,138
51,160,44
197,27,138
251,154,153
227,26,28
253,191,111
255,127,0
202,178,214
106,61,154
255,255,153
225,225,225
53,151,143
140,81,10]./255;

cond = {'aerobic','anaerobic'};
for i = 1:length(cond)
    figure
    hold on
for j = 1:length(sub)
idx = find(strcmpi(all_data(:,2),sub(j)) & strcmpi(all_data(:,4),cond{i}));
plot(dp(idx,1),dp(idx,2),'o','LineWidth',0.75,'MarkerSize',3,'Color',color(j,:));
end

set(gca,'FontSize',6,'FontName','Helvetica');
legend(sub,'Fontsize',6)
limit = max(dp(:));
xlim([0 limit]);
ylim([0 limit]);
x = [0:0.1:limit];
plot(x,x,'--k','LineWidth',0.75)
set(gca,'FontSize',6,'FontName','Helvetica');

ylabel('Simulated growth rate [1/h]','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('Experimental growth rate [1/h]','FontSize',7,'FontName','Helvetica','Color','k');

set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
'XAxisLocation','top', 'XTick', [],...
'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])
end
