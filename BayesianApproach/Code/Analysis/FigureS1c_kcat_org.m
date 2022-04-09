% to generate the kcat_orginsim figure for FigureS1

fid = fopen('Kcat_combination_wildtype_mutant_withUniprot.json', 'r');
data = textscan(fid,'%s','TreatAsEmpty','NA');
data = join(data{1,1}','');
allcombined = jsondecode(data{1});
fclose(fid);

alldata = [{allcombined.Organism}',{allcombined.ECNumber}',{allcombined.Value}',{allcombined.Type}',{allcombined.Substrate}'];


% type 
species = unique(alldata(:,1));
for i = 1:length(species)
    
idx = strcmp(alldata(:,4),'wildtype') & strcmp(alldata(:,1),species(i));
x(i,1) = sum(idx);
idx = strcmp(alldata(:,4),'mutant') & strcmp(alldata(:,1),species(i));
x(i,2) = sum(idx);
idx = strcmp(alldata(:,1),species(i));
x(i,3) = sum(idx);

end
[a,b] = sort(x(:,1),'descend');
species = species(b);
x = x(b,:);
bar([x(1:20,1:2);sum(x(21:end,1:2))],'stacked','LineWidth',0.5,'BarWidth',0.5);
xticks([1:1:21]);
xticklabels([species(1:20);{'Others'}])
xtickangle(90)
set(gcf,'position',[0 200 250 150]);
set(gca,'position',[0.2 0.5 0.6 0.4]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('{\itk}_c_a_t number','FontSize',7,'FontName','Helvetica','Color','k');
leg = legend({'Wildtype','Mutant'},'Location','northwest','Fontsize',6,'FontName','Helvetica');
leg.ItemTokenSize = [0,5];
legend('boxoff')