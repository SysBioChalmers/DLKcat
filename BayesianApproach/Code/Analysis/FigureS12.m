% RMSE new figure

% kcatCoverage Figure 5a b
% this function alculate the kcat coverage foe the enzymes in the model
% enzyme coverage
load('strains.mat')
current_path = pwd;
cd ../../Results/model_build_files/

i = 111
cd('model_classic')
z = load([strains{i},'_classic.mat']);
model = z.model;
enzymedata = z.enzymedata;
coverage(i,3) = length(find(~cellfun(@isempty,model.grRules)));
strain = z.strain;
growthdata = z.growthdata;
coverage(i,1) = length(enzymedata.kcat);

%if length(enzymedata.proteins) == length(model.genes)
res_enzyme(i,3) = length(model.genes);
%else
%     warning(['check enzymedata.proteins field num with model.genes field for species: ', species{i}]);
%end
a = join(enzymedata.enzyme, ' and ');
a = split(a,' and ');
allgene = unique(a); % enzyme that are included in the ecmodel
res_enzyme(i,1) = length(allgene);

cd ../model_dl
z = load([strains{i},'_dl.mat']);
model = z.model;
enzymedata = z.enzymedata;
coverage(i,2) = length(enzymedata.kcat);
a = join(enzymedata.enzyme, ' and ');
a = split(a,' and ');
allgene = unique(a);
res_enzyme(i,2) = length(allgene);
cd ../

coverage = coverage(:,1:2)./coverage(:,3);

figure
fig3 = bar([coverage(111,1),coverage(111,2)],0.5,'FaceColor',[33, 102, 172]/255,'EdgeColor',[33, 102, 172]/255,'LineWidth',0.75);
fig3.FaceAlpha = 0.3;
set(gca,'FontSize',6,'FontName','Helvetica');
xticklabels({'Mean value','Variance',7,'FontName','Helvetica','Color','k'})
ylabel(['{\itk}_c_a_t coverage for enzymes'],'FontSize',7,'FontName','Helvetica','Color','k');
xtickangle(30)
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
box off
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
    'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
    'XAxisLocation','top', 'XTick', [],...
    'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])

res_enzyme = res_enzyme(:,1:2)./res_enzyme(:,3);
figure
fig3 = bar([res_enzyme(111,1),res_enzyme(111,2)],0.5,'FaceColor',[33, 102, 172]/255,'EdgeColor',[33, 102, 172]/255,'LineWidth',0.75);
fig3.FaceAlpha = 0.3;
set(gca,'FontSize',6,'FontName','Helvetica');
xticklabels({'Mean value','Variance',7,'FontName','Helvetica','Color','k'})
ylabel(['{\itk}_c_a_t coverage for enzymes'],'FontSize',7,'FontName','Helvetica','Color','k');
xtickangle(30)
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
box off
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
    'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
    'XAxisLocation','top', 'XTick', [],...
    'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])


% c

cd ../../Results/model_build_files/
cd('model_classic')
i
z = load([strains{i},'_classic.mat']);
model = z.model;
enzymedata = z.enzymedata;
strain = z.strain;
growthdata = z.growthdata;
max_growth = z.max_growth;
rxn2block = z.rxn2block;
[~,tot_prot_weight,~,~,~,~,~,~] = sumBioMass(model);
tot_prot_weight = tot_prot_weight*0.5;

[rmse_final2,exp,simulated,growthdata,max_growth] = abc_matlab_max(model,enzymedata,enzymedata.kcat,tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);
result(i,1) = rmse_final2;

figure 
meanerror = sum((exp(:,1)-simulated(:,1)).^2)/length(exp(:,1));
exp_flux = exp(:,2:end);
exp_flux = exp_flux(:);
simulated_flux = simulated(:,2:end);
simulated_flux = simulated_flux(:);
simulated_flux(isnan(exp_flux)) = [];
exp_flux(isnan(exp_flux)) = [];

meanerror = sum((exp_flux(:,1)-simulated_flux(:,1)).^2)/length(exp_flux(:,1));
hold on
cd ../model_dl
z = load([strains{i},'_dl.mat']);
model = z.model;
enzymedata = z.enzymedata;

[rmse_final2,exp,simulated3,growthdata,max_growth] = abc_matlab_max(model,enzymedata,enzymedata.kcat,tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);
result(i,2) = rmse_final2;
cd ../model_bayesian
cd(strains{i})
a = dir('kcat_genra*.txt');
if ~isempty(a)
    a = {a.name};
    tmp = strrep(a,'kcat_genra','');
    tmp = strrep(tmp,'.txt','');
    tmp = cellfun(@str2double,tmp);
    tmp = max(tmp);
    a = ['kcat_genra',num2str(tmp),'.txt'];
    tmp = readmatrix(a,'FileType','text','Delimiter',',');
    tot_prot_weight = tmp(end-1,1);
    theta_100 = tmp(end,:); % is the rmse error
    result(i,3) = mean(theta_100);
    coverage(i,3) = coverage(i,2);
else
    result(i,3) = result(i,2);
    coverage(i,3) = coverage(i,2);
end
figure
fig3 = bar([result(111,1),result(111,2),result(111,3)],0.5,'FaceColor',[33, 102, 172]/255,'EdgeColor',[33, 102, 172]/255,'LineWidth',0.75);
fig3.FaceAlpha = 0.3;
set(gca,'FontSize',6,'FontName','Helvetica');
xticklabels({'Mean value','Variance',7,'FontName','Helvetica','Color','k'})
ylabel('RMSE','FontSize',7,'FontName','Helvetica','Color','k');
xtickangle(30)
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
box off
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
    'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
    'XAxisLocation','top', 'XTick', [],...
    'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])

cd ../../
