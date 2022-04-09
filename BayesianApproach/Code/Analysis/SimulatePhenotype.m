% SimulatePhenotype
%load species
current_path = pwd;
[~,~,growthrates] = xlsread('growthratedata.xlsx','growthrates'); % DLKcatBayesianApproach/Data
growthrates = growthrates(2:end,:);
species = unique(growthrates(:,4));
species = strrep(species,' ','_');

load('strains.mat');
species = intersect(species,strains);
cd ../../Results/model_build_files/
for i = 1:length(species)
    spec = species{i};
    disp([spec,':', num2str(i)])
    
    cd model_dl
    z = load([spec,'_dl.mat']);
    enzymedata = z.enzymedata;
    model = z.model;
    growthdata = z.growthdata;
    max_growth = z.max_growth;
    enzymedata.proteins = z.MWdata.genes;
    enzymedata.proteinMW = z.MWdata.MW;
    rxn2block = z.rxn2block;
    strain = z.strain;
    cd ../model_bayesian
    
    cd(spec)
    
    % plot rmse change
    nfound = length(dir('kcat_genra*.txt'));
    if nfound > 0
        tmp = readmatrix(['kcat_genra',num2str(nfound),'.txt'],'FileType','text','Delimiter',',');
        kcat_posterior = tmp(1:end-2,:);
        tot_prot_weight = tmp(end-1,1);
        ss = num2cell(kcat_posterior',1);
        [a,b] = arrayfun(@updateprior,ss);
        enzymedata.kcat = a';
    else
        [~,tot_prot_weight,~,~,~,~,~,~] = sumBioMass(model);
        tot_prot_weight = tot_prot_weight*0.5;
    end

    % save simulation
     [~,exp,simulated,growthdata,max_growth] = abc_matlab_max(model,enzymedata,enzymedata.kcat,tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);

    save([strain,'_sim_phen.mat'],'simulated','growthdata','max_growth')
    clearvars simulated_meadian simulated kcat_posterior
    cd ../../
end
cd(current_path)