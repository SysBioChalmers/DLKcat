function SimulatePhenotype(species)
for i = 1:length(species)
    spec = species{i};
    disp([spec,':', num2str(i)])
    
    cd ../KcatTuning/model_dl
    z = load([spec,'_dl.mat']);
    enzymedata = z.enzymedata;
    model = z.model;
    growthdata = z.growthdata;
    max_growth = z.max_growth;
    enzymedata.proteins = z.MWdata.genes;
    enzymedata.proteinMW = z.MWdata.MW;
    rxn2block = z.rxn2block;
    strain = z.strain;
    [~,tot_prot_weight,~,~,~,~,~,~] = sumBioMass(model);
    tot_prot_weight = tot_prot_weight*0.5;
    cd ../model_bayesian
    
    cd(spec)
    
    % plot rmse change
    try
        nfound = length(dir('kcat_genra*.txt'));
        tmp = readmatrix(['kcat_genra',num2str(nfound),'.txt'],'FileType','text','Delimiter',',');
        tot_prot_weight = tmp(end-1,1);
        theta_100 = tmp(end,:); % is the rmse error
        kcat_100 = tmp(1:end-2,:);
    catch
        warning('no posterior data found for species',strain)
    end
    % save simulation
    if exist('kcat_100','var')
        for m = 1:1
            disp([num2str(m),'/100'])
            [~,exp,simulated(:,:,m),growthdata,max_growth] = abc_matlab_max(model,enzymedata,kcat_100(:,m),tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);
        end
        simulated_meadian = median(simulated,3);
    else
        [~,exp,simulated(:,:,1),growthdata,max_growth] = abc_matlab_max(model,enzymedata,enzymedata.kcat,tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);
        simulated_meadian = median(simulated,3);
    end
    save([strain,'sim_phen.mat'],'simulated','growthdata','max_growth','simulated_meadian')
    clearvars simulated_meadian simulated kcat_100
    cd ../../
end