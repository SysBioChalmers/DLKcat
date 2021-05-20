

load('Results/kcat_pca_dl.mat')
pca_dl = pcaresult;
load('Results/kcat_pca_auto.mat')
pca_auto = pcaresult;   
rxnlist = intersect(pca_auto.rxnlist,pca_dl.rxnlist);

clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
[~,idx] = ismember(pca_dl.species,Strain_information(:,1));
pca_dl.clade = Strain_information(idx,2);

[~,idx] = ismember(rxnlist,pca_auto.rxnlist);
allpca.kcat = pca_auto.kcat(idx,:);
[~,idx] = ismember(rxnlist,pca_dl.rxnlist);
allpca.kcat = [allpca.kcat, pca_dl.kcat(idx,:)];
allpca.kcat (allpca.kcat == 0 ) = nan;
    all_norm = log10(allpca.kcat);
    all_norm(isnan(all_norm)) = 0;
  
    [~, score, ~, ~, explained, ~] = pca(all_norm','NumComponents',2);
  
   hold on;
for i = 1:13
idx = ismember(pca_dl.clade,clades(i));
h(i+1) = scatter(score(idx,1),score(idx,2),50,'o','filled','LineWidth',1);legend off;
end
legend(strrep(clades,'_',' \_'),'FontSize',5)
xlabel(['PC1 (',num2str(round(explained(1),1)),'%)'],'FontSize',7,'FontName','Helvetica','Color','k');
ylabel(['PC2 (',num2str(round(explained(2),1)),'%)'],'FontSize',7,'FontName','Helvetica','Color','k');
    
    
    