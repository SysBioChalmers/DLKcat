% this function is to identify crabtree positive species, do they have some
% significant kcat difference
inputpath = '/Users/feiranl/Documents/GitHub/MLKcat/ComplementaryScripts/KcatTuning/model_dl';
crabtree = {'Dekkera_bruxellensis';'Hanseniaspora_vinae';'Candida_glabrata';'Lachancea_fermentati';...
    'Lachancea_kluyveri';'Lachancea_thermotolerans';'Lachancea_waltii';'Nakaseomyces_castellii';'Nakaseomyces_delphensis';'Naumovozyma_castellii';...
    'Saccharomyces_cerevisiae';'Saccharomyces_eubayanus';'Saccharomyces_kudriavzevii';'Saccharomyces_mikatae';'Saccharomyces_paradoxus';...
    'Saccharomyces_uvarum';'Tetrapisispora_phaffii';'Torulaspora_delbrueckii';'Vanderwaltozyma_polyspora';'Zygosaccharomyces_bailii';'Zygosaccharomyces_rouxii';...
    'yHMPu5000034876_Tetrapisispora_iriomotensis';'yHMPu5000034710_Kluyveromyces_dobzhanskii';'yHMPu5000026152_Torulaspora_franciscae';'Schizosaccharomyces_pombe'};
nocrabtree = {'Kluyveromyces_marxianus','Lipomyces_starkeyi','yHMPu5000034761_Lipomyces_lipofer','yHMPu5000034760_Lipomyces_kononenkoae','yHMPu5000034758_Lipomyces_japonicus','yHMPu5000034757_Lipomyces_doorenjongii','yHMPu5000034754_Lipomyces_arxii','yHMPu5000034749_Lipomyces_mesembrius','yHMPu5000034748_Lipomyces_oligophaga','yHMPu5000034742_Lipomyces_suomiensis','Tortispora_caseinolytica','yHMPu5000035650_Trigonopsis_variabilis','yHMPu5000035282_Trigonopsis_vinaria','yHMPu5000034655_Botryozyma_nematodophila','Candida_infanticola','Saprochaete_clavata','Starmerella_bombicola_JCM9596','Wickerhamiella_domercqiae','Yarrowia_deformans','Yarrowia_keelungensis','Yarrowia_lipolytica','yHMPu5000035665_Middelhovenomyces_tepae','yHMPu5000035633_Candida_hispaniensis','yHMPu5000034674_Blastobotrys_muscicola','yHMPu5000034667_Blastobotrys_serpentis','yHMPu5000034661_Dipodascus_albidus','yHMPu5000034646_Wickerhamiella_cacticola','yHMPu5000034605_Spencermartinsiella_europaea','Alloascoidea_hylecoeti','yHMPu5000034604_Sporopachydermia_lactativora','Candida_arabinofermentans','Candida_boidinii_JCM9604','Komagataella_pastoris','Pichia_membranifaciens','yHMPu5000035675_Kregervanrija_fluxuum','yHMPu5000034904_Ogataea_nonfermentans','yHMPu5000034901_Ogataea_methylivora','yHMPu5000034893_Ogataea_philodendra','yHMPu5000034887_Ogataea_trehaloabstinens','yHMPu5000034637_Ogataea_populiabae','yHMPu5000034627_Pichia_heedii','yHMPu5000034625_Pichia_kudriavzevii','yHMPu5000026124_Ogataea_henricii','Babjeviella_inositovora','Candida_albicans','Candida_parapsilosis','Candida_tropicalis','Debaryomyces_hansenii','Priceomyces_haplophilus','Scheffersomyces_stipitis','Spathaspora_passalidarum','yHMPu5000041678_Debaryomyces_prosopidis','yHMPu5000035297_Priceomyces_castillae','yHMPu5000035296_Priceomyces_carsonii','yHMPu5000034999_Cephaloascus_fragrans','yHMPu5000034606_Priceomyces_medius','Ascoidea_rubescens','Wickerhamomyces_anomalus','yHMPu5000035673_Candida_orba','yHMPu5000035672_Phaffomyces_thermotolerans','yHMPu5000035671_Phaffomyces_antillensis','yHMPu5000035670_Phaffomyces_opuntiae','yHMPu5000035658_Starmera_amethionina','yHMPu5000035639_Wickerhamomyces_canadensis','yHMPu5000035286_Candida_azyma','yHMPu5000035274_Wickerhamomyces_alni','yHMPu5000035261_Candida_ponderosae','yHMPu5000035048_Barnettozyma_salicaria','yHMPu5000035046_Barnettozyma_populi','yHMPu5000035037_Candida_montana','Hanseniaspora_uvarum','Eremothecium_coryli','Eremothecium_cymbalariae','Eremothecium_gossypii','Eremothecium_sinecaudum','Kluyveromyces_lactis','yHMPu5000034709_Kluyveromyces_aestuarii'};
load('panmodel.mat')
model = panmodel;
% sol the model first by maxmize ATP yield
model = changeMedia(model,'D-glucose','MIN');
model_tmp = changeRxnBounds(model,'r_1714',-0.05,'l');
model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
model_tmp = changeObjective(model_tmp,'r_4046',1);
sol = optimizeCbModel(model_tmp,'max','one');
rxnTarget = model.rxns(find(sol.x));
model = panmodel;
model = changeMedia(model,'D-glucose','MIN');
model_tmp = changeRxnBounds(model,'r_1714',-0.05,'l');
model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
model_tmp = changeRxnBounds(model_tmp,'EX_protein_pool',-1000,'l');

model_tmp = changeObjective(model_tmp,'r_1761',1);
sol = optimizeCbModel(model_tmp,'max','one');

rxnTarget = [rxnTarget;model.rxns(find(sol.x))];

[siginificantEnzyme,pvalue,kcatresult,pathwayEnzyme] = kcatSignAnalysis(crabtree,nocrabtree,inputpath,0.01,rxnTarget);


% plot
rxnList = {'r_0300','r_0962'};
id = {'CIT','PYK'};
for i = 1:length(rxnList)
    group = [];
    idx = find(contains(kcatresult.rxns,rxnList(i)));
    group = [group;ones(length(crabtree),1);repmat(2,length(nocrabtree),1)];
    figure
    h = boxplot(kcatresult.kcat(idx,:)/3600',group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[197,27,138]/255);
    [r,p] = ttest2(kcatresult.kcat(idx,1:length(crabtree))./3600,kcatresult.kcat(idx,end-length(nocrabtree)+1:end)./3600);
    
    p
    text(1,300,['p = ' num2str(round(p,3))],'FontSize',6,'FontName','Helvetica','Color','k')
    ylabel(['Predicted {\itk}_c_a_t [1/s] for ',id{i}],'FontSize',7,'FontName','Helvetica','Color','k');
    set(gca,'position',[0.2 0.2 0.6 0.6]);
    set(gcf,'position',[0 200 150 150]);
    xticklabels({'Crabtree-positive','Crabtree-negative'})
    xtickangle(30)
    set(gca,'FontSize',6,'FontName','Helvetica');
    box off;
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
        'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
        'XAxisLocation','top', 'XTick', [],...
        'YAxisLocation','right', 'YTick', []);
    linkaxes([ax1, ax2])
    
end