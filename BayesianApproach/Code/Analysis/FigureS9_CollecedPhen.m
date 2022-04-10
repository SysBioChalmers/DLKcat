% This is to plot the figure S9

%% Figure Sa
fid2 = fopen('../../Data/343_phenotype_clade.tsv');
format = '%s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
    Strain_information(:,i) = temp{i};
end
fclose(fid2);

clades = unique(Strain_information(:,2));
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
group = [];
clade_av = [];

for i = 1:length(clades)
    idx = ismember(Strain_information(:,2),clades(i));
    clade_species(i,1) = length(find(idx)); % species number in each clade
end

% with data
load('species_withdata.mat')
[~,idx] = ismember(species_withdata,Strain_information(:,1));
x = Strain_information(idx,:);
for i = 1:length(clades)
idx = ismember(x(:,2),clades(i));
clade_species(i,2) = length(find(idx)); % species number in each clade
end

h = bar(clade_species,'stacked', 'FaceColor','flat','LineWidth',1);
h(1).CData = [178, 24, 43]/255;
h(2).CData = [33, 102, 172]/255;
h(1).FaceAlpha = 0.3;
h(2).FaceAlpha = 0.3;
h(1).BarWidth = 0.5;
h(2).BarWidth = 0.5;
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('No. species','FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'XTickLabel',[])
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
% remove the tick from right axis and top axis
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
           'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
           'XAxisLocation','top', 'XTick', [],... 
           'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])


%% Figure Sb
[~,~,growthrates] = xlsread('growthratedata.xlsx','growthrates');
growthrates = growthrates(2:end,:);
species_withdata = unique(growthrates(:,4));

a = tabulate(growthrates(:,7))
[~,idx] = sort(cell2mat(a(:,2)),'descend');
a_sorted = a(idx,:)

bar(cell2mat(a_sorted(:,2)),'LineWidth',0.5,'BarWidth',0.5,'FaceAlpha',0.3,'FaceColor',[178, 24, 43]/255);
xticks([1:1:16]);
xticklabels(a_sorted(:,1))
xtickangle(90)
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.3 0.6 0.6]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Entry number','FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('No. species','FontSize',7,'FontName','Helvetica','Color','k');
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
           'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
           'XAxisLocation','top', 'XTick', [],... 
           'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])

