function [rmse_final,exp,simulated,growthdata,max_growth]=abc_matlab_max(model,enzymedata,kcat_random_all,tot_prot_weight,growthdata,max_growth,proc,sample_generation,j,rxn2block)
nstep = sample_generation/proc;
rmse_final = zeros(1,nstep);
kcat_sample = kcat_random_all(:,(j-1)*nstep+1:j*nstep);


% get carbonnum for each exchange rxn to further calculation of error
if ~isfield(model,'excarbon')
    model = addCarbonNum(model);
end

for k = 1:nstep
    disp(['nstep:',num2str(k),'/',num2str(nstep)])
    kcat_random  = kcat_sample(:,k);
    %kcat_random = normrnd(enzymedata.kcat,enzymedata.kcat_var);
    prot_cost_info.value = enzymedata.MW./kcat_random;
    prot_cost_info.id = enzymedata.rxn_list;
    
    %% first search with substrate constrain
    objective = 'r_2111';
    osenseStr = 'max';
    if ~isempty(growthdata)
        [rmse_1,exp_1,simulated_1] = rmsecal(model,growthdata,true,objective,osenseStr,prot_cost_info,tot_prot_weight,rxn2block);
    else
        rmse_1 = [];
        exp_1 = [];
        simulated_1 = [];
    end
    %% second search for maxmial growth rate without constrain
    if ~isempty(max_growth)  % simulate the maximal growth rate
        [rmse_2,exp_2,simulated_2] = rmsecal(model,max_growth,false,objective,osenseStr,prot_cost_info,tot_prot_weight,rxn2block);
    else
        rmse_2 = [];
        exp_2 = [];
        simulated_2 = [];
    end
    exp = [exp_1;exp_2];
    simulated = [simulated_1;simulated_2];
    rmse_final(1,k) = mean([rmse_1,rmse_2],'omitnan');
    
    %% only output simulated result for one generation
    if nstep ~= 1 || sample_generation ~= 1
        simulated = [];
        exp = [];
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rmse,exp,simulated] = rmsecal(model,data,constrain,objective,osenseStr,prot_cost_info,tot_prot_weight,rxn2block)
exp_d = [];
sim = [];
simulated = zeros(length(data(:,1)),9);

for i = 1:length(data(:,1))
    exp = cell2mat(data(:,3:11)); % u sub ace eth gly pyr ethyl_acetate co2 o2
    exp = exp.*[1,-1,1,1,1,1,1,1,-1];
    ex_mets = {'growth',[data{i,2},' exchange'],'acetate exchange','ethanol exchange','glycerol exchange','pyruvate exchange','ethyl acetate exchange','carbon dioxide exchange','oxygen exchange'};
    [~,idx] = ismember(ex_mets,model.rxnNames);
    model_tmp = model;
    
    model_tmp = changeMedia(model_tmp,'D-glucose',data(i,16));
    model_tmp = changeRxnBounds(model_tmp,'r_1634',0,'b');
    model_tmp = changeRxnBounds(model_tmp,'r_1631',0,'b');

    if strcmp(data(i,14),'anaerobic') ||strcmp(data(i,14),'limited') 
        model_tmp = anaerobicModel(model_tmp);
    end
    if strcmp(data(i,14),'limited') 
         model_tmp.lb(strcmp(model_tmp.rxnNames,'oxygen exchange')) = -5;
    end
    if ~constrain
        model_tmp.lb(strcmp(model_tmp.rxns,'r_1714')) = 0;
        model_tmp.lb(strcmp(model_tmp.rxns,model.rxns(idx(2)))) = -1000; % not constrain the substrate usage
    else
        model_tmp.lb(strcmp(model_tmp.rxns,'r_1714')) = 0;
        model_tmp.lb(idx(2)) = exp(i,2);
    end
    
    sol_tmp = solveModel(model_tmp,objective,osenseStr,prot_cost_info,tot_prot_weight,'ibm_cplex');
    sol(:,i) = sol_tmp.x;
    
    tmp = ~isnan(exp(i,:));
    excarbon = model.excarbon(idx);
    excarbon(excarbon == 0) = 1;
    exp_tmp = exp(i,tmp).*excarbon(tmp);
    simulated_tmp = sol(idx(tmp),i)'.*excarbon(tmp); % normalize the growth rate issue by factor 10
    
    exp_block = zeros(1,length(setdiff(rxn2block,model_tmp.rxns(idx(2))))); % all zeros for blocked exchange mets exchange
    rxnblockidx = ismember(model_tmp.rxns,setdiff(rxn2block,model_tmp.rxns(idx(2))));
    simulated_block = sol(rxnblockidx,i)'.* model.excarbon(rxnblockidx); %
    exp_block = exp_block(simulated_block~=0);
    simulated_block = simulated_block(simulated_block~=0);
    if constrain
        rmse_tmp(i) = sqrt(immse([exp_tmp,exp_block], [simulated_tmp,simulated_block]));
    else
        if length(exp_tmp) >= 2
            rmse_tmp(i) = sqrt(immse(exp_tmp(1:2), simulated_tmp(1:2)));
        else
            rmse_tmp(i) = sqrt(immse(exp_tmp(1), simulated_tmp(1)));
        end
        
    end
    simulated(i,:) = sol(idx,i)';
end
rmse = sum(rmse_tmp)/length(data(:,1));
end
