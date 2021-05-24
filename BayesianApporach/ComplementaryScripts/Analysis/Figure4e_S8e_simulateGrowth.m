function Figure4e_S8e_simulateGrowth(species,carbonsource,aerobic)
if nargin < 3
    aerobic = 1;
end
if nargin < 2
    carbonsource = 'D-glucose';
end
current_path= pwd;
cd ../../Results/model_bayesian/

for i = 1:length(species)
    cd(species{i})
    % SimulateGrowth minmize glucose utake rate under different dilution rate
    dilutionrate = [0.05:0.05:0.9];
    
    %for m = 1:length(kcat_posterior(1,:))
    %for m = 1:100
       % load(['emodel_',species{i},num2str(m),'.mat']);
       m = 1
        load(['emodel_',species{i},'_Bayesian_DL_mean.mat']);
        if ~aerobic
            emodel = anaerobicModel(emodel);
        end
        for k = 1:length(dilutionrate)
            ex_mets = {'biomass pseudoreaction',[carbonsource,' exchange'],'acetate exchange','ethanol exchange','glycerol exchange','pyruvate exchange','ethyl acetate exchange','carbon dioxide exchange','oxygen exchange','EX_protein_pool'};
            [~,idx] = ismember(ex_mets,emodel.rxnNames);
            model_tmp = emodel;
            % first minimize carbon source uptake
            model_tmp.lb(strcmp(model_tmp.rxns,'r_1714')) = 0;
            model_tmp.lb(strcmp(model_tmp.rxns,model_tmp.rxns(idx(2)))) = -1000; % not constrain the substrate usage
            model_tmp.lb(strcmp(model_tmp.rxns,model_tmp.rxns(idx(1)))) = dilutionrate(k); % not constrain the substrate usage
            model_tmp = setParam(model_tmp,'obj',model_tmp.rxns(idx(2)),1);
            sol = optimizeCbModel(model_tmp);
            sol.f
            % fix the carbon source, minimize protein pool
            if ~isnan(sol.f) && sol.f ~= 0
                model_tmp = setParam(model_tmp,'lb',model_tmp.rxns(idx(2)),sol.f*1.00001);
                model_tmp = setParam(model_tmp,'lb',model_tmp.rxns(idx(10)),-1000);
                model_tmp = setParam(model_tmp,'obj',model_tmp.rxns(idx(10)),1);
                sol = optimizeCbModel(model_tmp);
                sol.f
                if ~isnan(sol.f) && sol.f ~= 0
                    sol_result(:,k,m) = sol.x;
                    simulated(k,:,m) = sol.x(idx)';
                end
            end
        end
        m
    end
    
    
    %plot exp using scatter
    cd ../../model_dl/
    z = load([species{i},'_dl.mat'])
   
    growthdata = z.growthdata;
    exp = growthdata(find(strcmpi(growthdata(:,2),carbonsource)& strcmpi(growthdata(:,14),'aerobic')),3:11)
    exp = abs(cell2mat(exp));
    color = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
     hold on
    for j = [1,3,7,8]
        scatter(exp(:,1),exp(:,j+1),30,'o','LineWidth',0.75,'MarkerFaceColor',color(j,:),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.3);
    end


% plot simulated
simulated = simulated.*[1,-1,1,1,1,1,1,1,-1,-0.001];
simulated_mean = median(simulated,3,'omitnan');
simulated_mean(all(simulated_mean==0,2),:) = [];
a = permute(simulated,[2 1 3]);
a = reshape(a,10,[])';
a(all(a==0,2),:) = [];
for i = 1:length(simulated_mean(:,1))
    idx = find(a(:,1) == simulated_mean(i,1))
    tmp = a(idx,:);
    tmp = sortrows(tmp)
    simulated_5th(i,:) = tmp(ceil(0.05*length(idx)),:);
    simulated_95th(i,:) = tmp(ceil(0.95*length(idx)),:);
end
hold on
k = 1
for j = [1,3,7,8]
    h(k) = plot(simulated_mean(:,1),simulated_mean(:,j+1),'-','LineWidth',0.75,'color',color(j,:));
    k = k+1;
    %     X_plot = [simulated_mean(:,1)', fliplr(simulated_mean(:,1)')];
    %     Y_plot = [simulated_5th(:,j+1)', fliplr(simulated_95th(:,j+1)')];
    %     %nanx = isnan(X_plot)| isnan(Y_plot)
    %     %X_plot = X_plot(~nanx)
    %     %Y_plot = Y_plot(~nanx)
    %     Y_plot(isnan(Y_plot)) = 0;
    %     fill(X_plot, Y_plot , 1,....
    %         'facecolor',color(j,:), ...
    %         'edgecolor','none', ...
    %         'facealpha', 0.3);
end
hold off
set(gca,'FontSize',6,'FontName','Helvetica');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
k = 0
for j = [1,3,7,8]
    text(0.05,7+k,ex_mets(j+1),'FontSize',6,'FontName','Helvetica','Color',color(j,:))
    k = k+0.8;
end
box on
ylabel('Exchange rate[mmol/gDW/h]','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('Dilution rate[1/h]')
xlim([0 0.4]);
saveas(h,['growthfigure_',species{i},'.pdf']);
%% max growth
for m = 1:length(kcat_posterior(1,:))
    load(['emodel_',strain,num2str(m),'.mat'])
    %emodel = changeGAM(emodel,24.6807,1.3972);
    %emodel = blockbyproduct(emodel);
    if ~isempty(max_growth)  % simulate the maximal growth rate
        exp = cell2mat(max_growth(:,3:10)); % u sub ace eth gly pyr co2 o2
        exp = exp.*[1,-1,1,1,1,1,1,-1];
        for j = 1:length(max_growth(:,1))
            ex_mets = {'biomass pseudoreaction',[max_growth{j,2},' exchange'],'acetate exchange','ethanol exchange','glycerol exchange','pyruvate exchange','ethyl acetate exchange','carbon dioxide exchange','oxygen exchange'};
            [~,idx] = ismember(ex_mets,emodel.rxnNames);
            model_tmp = emodel;
            if strcmp(max_growth(j,14),'anaerobic')
                model_tmp = anaerobicModel(model_tmp);
            end
            model_tmp.lb(strcmp(model_tmp.rxns,'r_1714')) = 0;
            model_tmp.lb(strcmp(model_tmp.rxns,model.rxns(idx(2)))) = -1000; % not constrain the substrate usage
            sol = optimizeCbModel(model_tmp);
            max_growth_simulated(j,m) = sol.f;
            max_sol_result(:,j,m) = sol.x;
        end
    end
    m
end
max_growth_ave = mean(max_growth_simulated,2);
max_growth_std = std(max_growth_simulated,0,2);
max_growth_ave = [cell2mat(max_growth(:,3)),max_growth_ave];
max_growth_std = [zeros(length(max_growth(:,3)),1),max_growth_std];

figure
z2 = bar(max_growth_ave);
hold on
% Calculate the number of bars in each group
nbars = size(max_growth_ave, 2);
% Get the x coordinate of the bars
x = [];
for m = 1:nbars
    x = [x ; z2(m).XEndPoints];
end
% Plot the errorbars
errorbar(x',max_growth_ave,max_growth_std,'k','Marker','none','LineStyle','none','LineWidth',0.3,'CapSize',4)
hold off

set(gca,'FontSize',8,'FontName','Helvetica');
set(gca,'ycolor','k');
color = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
for j = 1:nbars
    z2(j).EdgeColor = 'k';
    z2(j).FaceColor = color(j,:);
    z2(j).FaceAlpha = 0.2;
end
xticklabels(max_growth(:,2))
xtb = get(gca,'XTickLabel');
xt = get(gca,'XTick');
yt = get(gca,'YTick');
xtextp = xt;
ytextp = yt(1)*ones(1,length(xt));
text(xtextp,ytextp-0.03,xtb,'HorizontalAlignment','center','FontSize',8,'FontName','Helvetica');
set(gca,'xticklabel','');
ylabel('Maximal Growth rates[1/h]','FontSize',8,'FontName','Helvetica','Color','k');
legend({'exp','model'},'FontSize',5)
set(gcf,'position',[600 0 450 170]);
set(gca,'position',[0.1 0.4 0.27 0.56]);
saveas(z2,['growthfigure_',species{i},'.pdf']);
end
