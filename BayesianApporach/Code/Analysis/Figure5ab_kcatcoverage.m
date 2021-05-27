% Figure for kcat coverage comparsion

load('res_kcatcoveragedl.mat') % generate from the kcatCoverage.m
dl_res_rxn = res_rxn;
dl_res_enzyme = res_enzyme;
dl_final_rxn = dl_res_rxn(:,1)./dl_res_rxn(:,2);
dl_final_enzyme = dl_res_enzyme(:,1)./dl_res_enzyme(:,2);
dl_species = species;

load('res_kcatcoverageauto.mat')
auto_res_rxn = res_rxn;
auto_res_enzyme = res_enzyme;
auto_final_rxn = auto_res_rxn(:,1)./auto_res_rxn(:,2);
auto_final_enzyme = auto_res_enzyme(:,1)./auto_res_enzyme(:,2);
auto_species = species;

tmp = cellfun(@strcmp,dl_species,auto_species);
if all(tmp)
    figure
    hold on
    fig1(1) = cdfplot(auto_final_rxn);
    fig1(2)  = cdfplot(dl_final_rxn);
    ylabel('Ratio of ecModel Enzyme coverage','FontSize',8,'FontName','Helvetica','Color','k');
    xlabel('species number','FontSize',8,'FontName','Helvetica','Color','k');
    legend({'Auto','DL'},'FontSize',6)
    set(gcf,'position',[500 200 180 130]);
    set(gca,'position',[0.27 0.33 0.7 0.63]);
    saveas(fig1,'enzymecoverage.pdf');
    
    figure
     hold on
    fig2(1) = cdfplot(auto_final_enzyme);
    fig2(2)  = cdfplot(dl_final_enzyme);
    ylabel('Ratio of ecModel Enzyme coverage','FontSize',8,'FontName','Helvetica','Color','k');
    xlabel('species number','FontSize',8,'FontName','Helvetica','Color','k');
    legend({'Auto','DL'},'FontSize',6)
    set(gcf,'position',[500 200 180 130]);
    set(gca,'position',[0.27 0.33 0.7 0.63]);
    saveas(fig2,'rxncoverage.pdf');
end
