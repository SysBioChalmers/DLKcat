% This function is to calculate the protein efficieny for different yeast
% species

% load yeast species models
crabtree = {'Dekkera_bruxellensis';'Hanseniaspora_vinae';'Candida_glabrata';'Lachancea_fermentati';...
'Lachancea_kluyveri';'Lachancea_thermotolerans';'Lachancea_waltii';'Nakaseomyces_castellii';'Nakaseomyces_delphensis';'Naumovozyma_castellii';...
'Saccharomyces_cerevisiae';'Saccharomyces_eubayanus';'Saccharomyces_kudriavzevii';'Saccharomyces_mikatae';'Saccharomyces_paradoxus';...
'Saccharomyces_uvarum';'Tetrapisispora_phaffii';'Torulaspora_delbrueckii';'Vanderwaltozyma_polyspora';'Zygosaccharomyces_bailii';'Zygosaccharomyces_rouxii';...
'yHMPu5000034876_Tetrapisispora_iriomotensis';'yHMPu5000034710_Kluyveromyces_dobzhanskii';'yHMPu5000026152_Torulaspora_franciscae';'Schizosaccharomyces_pombe'};
nocrabtree = {'Kluyveromyces_marxianus','Lipomyces_starkeyi','yHMPu5000034761_Lipomyces_lipofer','yHMPu5000034760_Lipomyces_kononenkoae','yHMPu5000034758_Lipomyces_japonicus','yHMPu5000034757_Lipomyces_doorenjongii','yHMPu5000034754_Lipomyces_arxii','yHMPu5000034749_Lipomyces_mesembrius','yHMPu5000034748_Lipomyces_oligophaga','yHMPu5000034742_Lipomyces_suomiensis','Tortispora_caseinolytica','yHMPu5000035650_Trigonopsis_variabilis','yHMPu5000035282_Trigonopsis_vinaria','yHMPu5000034655_Botryozyma_nematodophila','Candida_infanticola','Saprochaete_clavata','Starmerella_bombicola_JCM9596','Wickerhamiella_domercqiae','Yarrowia_deformans','Yarrowia_keelungensis','Yarrowia_lipolytica','yHMPu5000035665_Middelhovenomyces_tepae','yHMPu5000035633_Candida_hispaniensis','yHMPu5000034674_Blastobotrys_muscicola','yHMPu5000034667_Blastobotrys_serpentis','yHMPu5000034661_Dipodascus_albidus','yHMPu5000034646_Wickerhamiella_cacticola','yHMPu5000034605_Spencermartinsiella_europaea','Alloascoidea_hylecoeti','yHMPu5000034604_Sporopachydermia_lactativora','Candida_arabinofermentans','Candida_boidinii_JCM9604','Komagataella_pastoris','Pichia_membranifaciens','yHMPu5000035675_Kregervanrija_fluxuum','yHMPu5000034904_Ogataea_nonfermentans','yHMPu5000034901_Ogataea_methylivora','yHMPu5000034893_Ogataea_philodendra','yHMPu5000034887_Ogataea_trehaloabstinens','yHMPu5000034637_Ogataea_populiabae','yHMPu5000034627_Pichia_heedii','yHMPu5000034625_Pichia_kudriavzevii','yHMPu5000026124_Ogataea_henricii','Babjeviella_inositovora','Candida_albicans','Candida_parapsilosis','Candida_tropicalis','Debaryomyces_hansenii','Priceomyces_haplophilus','Scheffersomyces_stipitis','Spathaspora_passalidarum','yHMPu5000041678_Debaryomyces_prosopidis','yHMPu5000035297_Priceomyces_castillae','yHMPu5000035296_Priceomyces_carsonii','yHMPu5000034999_Cephaloascus_fragrans','yHMPu5000034606_Priceomyces_medius','Ascoidea_rubescens','Wickerhamomyces_anomalus','yHMPu5000035673_Candida_orba','yHMPu5000035672_Phaffomyces_thermotolerans','yHMPu5000035671_Phaffomyces_antillensis','yHMPu5000035670_Phaffomyces_opuntiae','yHMPu5000035658_Starmera_amethionina','yHMPu5000035639_Wickerhamomyces_canadensis','yHMPu5000035286_Candida_azyma','yHMPu5000035274_Wickerhamomyces_alni','yHMPu5000035261_Candida_ponderosae','yHMPu5000035048_Barnettozyma_salicaria','yHMPu5000035046_Barnettozyma_populi','yHMPu5000035037_Candida_montana','Hanseniaspora_uvarum','Eremothecium_coryli','Eremothecium_cymbalariae','Eremothecium_gossypii','Eremothecium_sinecaudum','Kluyveromyces_lactis','yHMPu5000034709_Kluyveromyces_aestuarii'};
species = [crabtree;nocrabtree'];
cd ../../Results/model_bayesian
for i = 1:length(species)
    cd(species{i})
    i
    for k = 1:1
        disp([num2str(i),' for species: ',species{i}])
        %load(['emodel_',species{i},num2str(k),'.mat']);
        load(['emodel_',species{i},'_Bayesian_DL_mean.mat']);
        %load(['emodel_',species{i},'_dl.mat']);
        model = emodel;
        % sol the model first by maxmize ATP yield
        model_tmp = changeRxnBounds(model,'r_1714',-0.05,'l');
        model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
        model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
        model_tmp = changeRxnBounds(model_tmp,'EX_protein_pool',-1000,'l');
        model_tmp = changeObjective(model_tmp,'r_4046',1);
        sol = optimizeCbModel(model_tmp);
        % fix the ATP, minimize protein pool
        if ~isnan(sol.f) && sol.f ~= 0
            HY_atp(k,i) = sol.f;
            model_tmp = changeRxnBounds(model_tmp,'r_4046',sol.f*1.0001,'u');
            model_tmp = changeRxnBounds(model_tmp,'r_4046',sol.f*0.9999,'l');
            model_tmp = changeObjective(model_tmp,'EX_protein_pool',1);
            sol = optimizeCbModel(model_tmp,'max');
            HY_protein(k,i) = sol.x(strcmp(model_tmp.rxns,'EX_protein_pool'));
        end
        % sol the model LY pathway by maxmize ethanol
        model_tmp = changeRxnBounds(model,'r_1714',-0.05,'l');
        model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
        model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
        model_tmp = changeRxnBounds(model_tmp,'EX_protein_pool',-1000000,'l');
     
        model_tmp = changeObjective(model_tmp,'r_1761',1);
        sol = optimizeCbModel(model_tmp);
        % fix the ATP, minimize protein pool
        if ~isnan(sol.f) && sol.f ~= 0
            model_tmp = changeRxnBounds(model_tmp,'r_1761',sol.f*1.0001,'u');
            model_tmp = changeRxnBounds(model_tmp,'r_1761',sol.f*0.9998,'l');
            model_tmp = changeObjective(model_tmp,'EX_protein_pool',1);
            sol = optimizeCbModel(model_tmp,'max','one');
            LY_protein(k,i) = sol.x(strcmp(model_tmp.rxns,'EX_protein_pool'));
            model_tmp = changeObjective(model_tmp,'r_4046',1);
            sol = optimizeCbModel(model_tmp,'max');
            LY_atp(k,i) = sol.f;
        end
    end
    cd ../
end
HY_protein_mean = abs(mean(HY_protein,1)./0.05);
HY_atp_mean = abs(mean(HY_atp,1)./0.05);
LY_atp_mean = abs(mean(LY_atp,1)./0.05);
LY_protein_mean = abs(mean(LY_protein,1)./0.05);

ratio_HY(:,1) = HY_atp_mean./HY_protein_mean;
ratio_LY(:,1) = LY_atp_mean./LY_protein_mean;
ratio(:,1) = ratio_LY./ratio_HY;

crabtree = {'Dekkera_bruxellensis';'Hanseniaspora_vinae';'Candida_glabrata';'Lachancea_fermentati';...
'Lachancea_kluyveri';'Lachancea_thermotolerans';'Lachancea_waltii';'Nakaseomyces_castellii';'Nakaseomyces_delphensis';'Naumovozyma_castellii';...
'Saccharomyces_cerevisiae';'Saccharomyces_eubayanus';'Saccharomyces_kudriavzevii';'Saccharomyces_mikatae';'Saccharomyces_paradoxus';...
'Saccharomyces_uvarum';'Tetrapisispora_phaffii';'Torulaspora_delbrueckii';'Vanderwaltozyma_polyspora';'Zygosaccharomyces_bailii';'Zygosaccharomyces_rouxii';...
'yHMPu5000034876_Tetrapisispora_iriomotensis';'yHMPu5000034710_Kluyveromyces_dobzhanskii';'yHMPu5000026152_Torulaspora_franciscae';'Schizosaccharomyces_pombe'};
nocrabtree = {'Kluyveromyces_marxianus';'Lipomyces_starkeyi';'yHMPu5000034761_Lipomyces_lipofer';'yHMPu5000034760_Lipomyces_kononenkoae';'yHMPu5000034758_Lipomyces_japonicus';'yHMPu5000034757_Lipomyces_doorenjongii';'yHMPu5000034754_Lipomyces_arxii';'yHMPu5000034749_Lipomyces_mesembrius';'yHMPu5000034748_Lipomyces_oligophaga';'yHMPu5000034742_Lipomyces_suomiensis';'Tortispora_caseinolytica';'yHMPu5000035650_Trigonopsis_variabilis';'yHMPu5000035282_Trigonopsis_vinaria';'yHMPu5000034655_Botryozyma_nematodophila';'Candida_infanticola';'Saprochaete_clavata';'Starmerella_bombicola_JCM9596';'Wickerhamiella_domercqiae';'Yarrowia_deformans';'Yarrowia_keelungensis';'Yarrowia_lipolytica';'yHMPu5000035665_Middelhovenomyces_tepae';'yHMPu5000035633_Candida_hispaniensis';'yHMPu5000034674_Blastobotrys_muscicola';'yHMPu5000034667_Blastobotrys_serpentis';'yHMPu5000034661_Dipodascus_albidus';'yHMPu5000034646_Wickerhamiella_cacticola';'yHMPu5000034605_Spencermartinsiella_europaea';'Alloascoidea_hylecoeti';'yHMPu5000034604_Sporopachydermia_lactativora';'Candida_arabinofermentans';'Candida_boidinii_JCM9604';'Komagataella_pastoris';'Pichia_membranifaciens';'yHMPu5000035675_Kregervanrija_fluxuum';'yHMPu5000034904_Ogataea_nonfermentans';'yHMPu5000034901_Ogataea_methylivora';'yHMPu5000034893_Ogataea_philodendra';'yHMPu5000034887_Ogataea_trehaloabstinens';'yHMPu5000034637_Ogataea_populiabae';'yHMPu5000034627_Pichia_heedii';'yHMPu5000034625_Pichia_kudriavzevii';'yHMPu5000026124_Ogataea_henricii';'Babjeviella_inositovora';'Candida_albicans';'Candida_parapsilosis';'Candida_tropicalis';'Debaryomyces_hansenii';'Priceomyces_haplophilus';'Scheffersomyces_stipitis';'Spathaspora_passalidarum';'yHMPu5000041678_Debaryomyces_prosopidis';'yHMPu5000035297_Priceomyces_castillae';'yHMPu5000035296_Priceomyces_carsonii';'yHMPu5000034999_Cephaloascus_fragrans';'yHMPu5000034606_Priceomyces_medius';'Ascoidea_rubescens';'Wickerhamomyces_anomalus';'yHMPu5000035673_Candida_orba';'yHMPu5000035672_Phaffomyces_thermotolerans';'yHMPu5000035671_Phaffomyces_antillensis';'yHMPu5000035670_Phaffomyces_opuntiae';'yHMPu5000035658_Starmera_amethionina';'yHMPu5000035639_Wickerhamomyces_canadensis';'yHMPu5000035286_Candida_azyma';'yHMPu5000035274_Wickerhamomyces_alni';'yHMPu5000035261_Candida_ponderosae';'yHMPu5000035048_Barnettozyma_salicaria';'yHMPu5000035046_Barnettozyma_populi';'yHMPu5000035037_Candida_montana';'Hanseniaspora_uvarum';'Eremothecium_coryli';'Eremothecium_cymbalariae';'Eremothecium_gossypii';'Eremothecium_sinecaudum';'Kluyveromyces_lactis';'yHMPu5000034709_Kluyveromyces_aestuarii'};
[~,idx] = ismember([crabtree;nocrabtree],species);
ratio_final = num2cell([ratio(idx),HY_atp_mean(idx)',HY_protein_mean(idx)',LY_atp_mean(idx)',LY_protein_mean(idx)',ratio_HY(idx),ratio_LY(idx)]);
ratio_id = species(idx);
ratio_final(:,8) = [repmat({'Crabtree'},length(crabtree),1);repmat({'Crabtree-negative'},length(nocrabtree),1)];
ratio_final(:,9) = ratio_id;

writecell(ratio_final,'../../Results/res_Crabtree.txt','Delimiter',',','QuoteStrings',false)

violin = violinplot(cell2mat(ratio_final(:,6))./cell2mat(ratio_final(:,7)),ratio_final(:,8),'ShowNotches',false,'ShowMean' ,false,'ViolinAlpha',1,'EdgeColor',[0,0,0],'ShowData',false,'BoxColor',[1,1,1]);
violin(1).ViolinColor = [253,224,221]./255;
violin(2).ViolinColor = [250,159,181]./255;
violin(2).MedianPlot.SizeData = 1;
violin(1).MedianPlot.SizeData = 1;
xticklabels({'Crabtree-positive','Crabtree-negative'})
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel(['Protein efficiency'],'FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
% remove the tick from right axis and top axis
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
           'Color','None','XColor','k','YColor','k', 'LineWidth', 1,...
           'XAxisLocation','top', 'XTick', [],... 
           'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])
