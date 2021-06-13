% Figure 5d
% This function is to plot the figure5d, which is mainly focusing on
% evalulation of the kcat tuning process
% Running this function requires that all tuned emodels are genrated

[~,~,growthrates] = xlsread('growthratedata.xlsx','growthrates');
growthrates = growthrates(2:end,:);
species = unique(growthrates(:,4));
species = strrep(species,' ','_');


% Figure 3A this plot is to show the growth rate prediction vs experimental
% values
result_constrain = [];
result_max = [];
cd ../../Results/model_Bayesian_max
for i = 1:length(species)
    cd(species{i})% later change to cd(species{i}) use species id to name the folder
    load([species{i},'sim_phen.mat']) % is generated in the function simulategrowth.m
    result_constrain = [result_constrain;growthdata(:,[1:3,14]),num2cell(simulated_meadian(1:length(growthdata(:,1)),1))];
    result_max = [result_max;max_growth(:,[1:3,14]),num2cell(simulated_meadian(length(growthdata(:,1))+1:end,1))];
    cd ../
end

% datapoint
all_data = [result_constrain;result_max];
dp = cell2mat(all_data(:,[3,5]));

sub = unique(all_data(:,2));
color = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153;2,2,2]./255;

hold on
for j = 1:length(sub)
idx = find(strcmpi(all_data(:,2),sub(j)) & strcmpi(all_data(:,4),'aerobic'));
plot(dp(idx,1),dp(idx,2),'o','LineWidth',0.75,'MarkerSize',5);
end
for j = 1:length(sub)
idx = find(strcmpi(all_data(:,2),sub(j)) & strcmpi(all_data(:,4),'anaerobic'));
plot(dp(idx,1),dp(idx,2),'<','LineWidth',0.75,'MarkerSize',5);
end
set(gca,'FontSize',6,'FontName','Helvetica');
legend(sub,'Fontsize',6)
limit = max(dp(:));
xlim([0 limit]);
ylim([0 limit]);
x = [0:0.1:limit];
plot(x,x,'--k','LineWidth',0.75)

ylabel('Simulated growth rate [1/h]','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('Experimental growth rate [1/h]','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[500 200 180 130]);
set(gca,'position',[0.27 0.33 0.7 0.63]);

