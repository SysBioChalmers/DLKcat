function TestBayesianModelGeneration
initcluster


species = {'Saccharomyces_cerevisiae'};

generation = 150;

cd('../../Results')
for i = 1:1
    cd('model_dl')
    mkdir([species{i},'_test'])
    z = load([species{i},'_dl.mat']);
    cd(species{i})
    enzymedata = z.enzymedata;
    max_growth = z.max_growth;
    growthdata = z.growthdata;
    model = z.model;
    strain = z.strain;
    growthdata_train = growthdata([2:2:length(growthdata(:,1))],:);
    growthdata_test = growthdata([1:2:length(growthdata(:,1))],:);
    max_growth_train = max_growth([2:2:length(max_growth(:,1))],:);
    max_growth_test = max_growth([1:2:length(max_growth(:,1))],:);
    %load('rxn2block.mat')
    rxn2block = z.rxn2block;
    [~,tot_prot_weight,~,~,~,~,~,~] = sumBioMass(model);
    tot_prot_weight = tot_prot_weight*0.5; % metabolic enzyme takes the 50%
    if strcmp(strain,'Kluyveromyces_marxianus')
        tot_prot_weight = 0.325; % from bionumber ID 117044
    elseif strcmp(strain,'Kluyveromyces_lactis')
        tot_prot_weight = 0.245; % from the published GEM k.lactis
    end
    
    proc = 18;
    numPerGeneration = 126;
    generation = 100;
    rejectnum = 0.2;
    
    if isempty(max_growth) && isempty(growthdata)
        max_growth = {strain,'D-glucose',0.2,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,'aerobic','Batch','MIN'};
        [~,~,simulated3,~,~] = abc_matlab_max(model,enzymedata,enzymedata.kcat,tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);
        if simulated3(1,1) > 0.2
            max_growth = {strain,'D-glucose',simulated3(1,1),nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,'aerobic','Batch','MIN'};
        end
    end
    
    D = abc_matlab_max(model,enzymedata,enzymedata.kcat,tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);
    D_100 = D;
    theta_100 = [];
    theta_100_test = [];
    kcat_100 = [];
    sampledgeneration = 1;
    while D > rejectnum && D_100 > 0.5 % range for D_1 to D_100 is within 0.3 and D_100 should < 0.5
        if continued
            nfound = length(dir('kcat_genra*.txt'));
            if nfound > 0 
            tmp = readmatrix(['kcat_genra',num2str(nfound),'.txt'],'FileType','text','Delimiter',',');
            theta_100 = tmp(end-1,:);
            theta_100_test = tmp(end,:);
            kcat_100 = tmp(1:end-3,:);
            tot_prot_weight = tmp(end-2,1);
            sampledgeneration = nfound + 1;
            % recalculate the sigma and mu
            ss = num2cell(kcat_100',1);
            [a,b] = arrayfun(@updateprior,ss);
            enzymedata.kcat = a';
            enzymedata.kcat_var = b';
            end
            continued = false; % close the loop after the first one
        end
        
        if sampledgeneration <= generation
            disp(['No.',num2str(sampledgeneration),' generation'])
            % generate a
            old = theta_100;
            old_test = theta_100_test;
            kcat_old_100 = kcat_100;
            %repeat a generation
            
            if sampledgeneration == 1
                sample_generation = 144;
            else
                sample_generation = numPerGeneration;
            end
            % generate one generation sample of kcats
            kcat_random_all = arrayfun(@getrSample, enzymedata.kcat,enzymedata.kcat_var,enzymedata.enzyme_ec_kcat_range(:,1),enzymedata.enzyme_ec_kcat_range(:,2),repmat(sample_generation,length(enzymedata.kcat),1),'UniformOutput',false);
            kcat_random_all = cell2mat(kcat_random_all);
            disp('kcat random finish')
            %start sampleing
            parfor j = 1:proc
                j
                rmse_final = abc_matlab_max(model,enzymedata,kcat_random_all,tot_prot_weight,growthdata_train,max_growth_train,proc,sample_generation,j,rxn2block)
                rmse_final_test = abc_matlab_max(model,enzymedata,kcat_random_all,tot_prot_weight,growthdata_test,max_growth_test,proc,sample_generation,j,rxn2block)
                new_tmp{j} = rmse_final;
                new_tmp_test{j} = rmse_final_test;
            end
            disp('rmse calculation finish')
            %find the best 100 samples
            %         for k = 1:length(new_tmp)
            %             tmp = new_tmp{k};
            %             new(:,(k-1)*proc+1:k*proc) = tmp;
            %         end
            new = cell2mat(new_tmp);
            new_test = cell2mat(new_tmp_test);
            theta = [new,old];
            theta_test = [new_test,old_test];
            kcat = [kcat_random_all,kcat_old_100];
            
            % Initialize an empty set to store the best 100  after each step
            [~,D_idx]= sort(theta,'ascend');
            theta_100 = theta(D_idx(1:100));
            theta_100_test = theta_test(D_idx(1:100));
            D = abs(theta_100(100)-theta_100(1)); % the largest one theta - smallest theta
            D_100 = theta_100(100);
            kcat_100 = kcat(:,D_idx(1:100));
            writematrix([kcat_100;repmat(tot_prot_weight,1,length(theta_100));theta_100;theta_100_test],['kcat_genra',num2str(sampledgeneration),'.txt'])
            
            % recalculate the sigma and mu
            ss = num2cell(kcat_100',1);
            [a,b] = arrayfun(@updateprior,ss);
            enzymedata.kcat = a';
            enzymedata.kcat_var = b';
            %enzymedata.kcat_var(enzymedata.kcat_var < 1) = 1;
            sampledgeneration = sampledgeneration + 1;
        else
            D = rejectnum;
            D_100 = D;
        end
    end
    cd('../')
end