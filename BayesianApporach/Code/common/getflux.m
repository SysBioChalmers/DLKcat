function [model,sol_result_mean] = getflux(strain,chemostat,media,anox,carbon,inputpath,state)

if ~strcmp(chemostat,'max')
    dilutionrate = str2double(chemostat);
    chemostat = 'chemo';
elseif strcmp(chemostat,'max')
    chemostat = 'batch';
end
current_path = pwd;
cd(inputpath)

if strcmp(state,'DL') || strcmp(state,'classic') || strcmp(state,'global')|| strcmp(state,'Posterior_mean')
    num = 1;
elseif strcmp(state,'Posterior')
    num = 5;
end
for m = 1:num
    m
    if strcmp(state,'DL')
        load(['emodel_',strain,'_DL.mat'])
    elseif strcmp(state,'Posterior')
        load(['emodel_',strain,num2str(m),'.mat'])
    elseif strcmp(state,'Posterior_mean')
        load(['emodel_',strain,'_Posterior_mean.mat'])
    elseif strcmp(state,'classic')
        load(['emodel_',strain,'_classic.mat'])
    elseif strcmp(state,'global')
        cd('../../model_dl')
        z = load([strain,'_dl.mat']);
        if strcmp(strain,'Kluyveromyces_marxianus')
            tot_prot = 0.325;
        elseif strcmp(strain,'Kluyveromyces_lactis')
             tot_prot = 0.245;
        else
             tot_prot = 0.23;
        end
        z.enzymedata.kcat(1:length(z.enzymedata.kcat)) = 79*3600;% median kcat of the metabolic enzyme Median of Central-CE (PMID: 21506553)
        emodel = convertToGeckoModel(z.model,z.enzymedata,tot_prot);
    end
    model = emodel;
    if strcmp(anox,'anaerobic')
        anox = 1;
    else
        anox = 0;
    end
    [model,~] = changeMedia(model,carbon,media,anox);
    if strcmpi(chemostat,'batch')
            [~,idx] = ismember(strcat(carbon,' exchange'),model.rxnNames);
            model = changeRxnBounds(model,'r_1714',0,'l');
            model = changeRxnBounds(model,model.rxns(idx),-1000,'l');
            model = changeObjective(model,model.rxns(strcmp(model.rxnNames,'growth')),1);
            try
            sol = optimizeCbModel(model,'max','one');
            if isempty(sol.f) | isnan(sol.f)
                try
                sol = optimizeCbModel(model);
                disp('loop allowed')
                catch
                    warning('no optimal solution for allowing loop')
                    sol.x = zeros(length(model.rxns),1);
                    sol.f = nan;
                end
            end
            catch
                disp('skip this one')
                sol.x = zeros(length(model.rxns),1);
                sol.f = nan;
            end
    else
            [~,idx] = ismember(strcat(carbon,' exchange'),model.rxnNames);
            model = changeRxnBounds(model,'r_1714',0,'l');
            model = changeRxnBounds(model,model.rxns(idx),-1000,'l');
            model = changeRxnBounds(model,'EX_protein_pool',-1000,'l'); % unlimited protein pool
            model = changeRxnBounds(model,model.rxns(strcmp(model.rxnNames,'growth')|strcmp(model.rxnNames,'biomass exchange')),dilutionrate);
            model = changeObjective(model,model.rxns(idx),1);
            try
            sol = optimizeCbModel(model,'max','one');
            catch
                disp('no optimal solution for without loop')
                sol.x = zeros(length(model.rxns),1);
                sol.f = nan;
            end
            if isempty(sol.f) | isnan(sol.f)
                try
                sol = optimizeCbModel(model);
                disp('loop allowed')
                catch
                    warning('no optimal solution for allowing loop')
                    sol.x = zeros(length(model.rxns),1);
                    sol.f = nan;
                end
            end
    end
   if isempty(sol.f) | isnan(sol.f)
     warning('no optimal solution for this case')
     sol.x = zeros(length(model.rxns),1);
   end

    sol_result(:,m) = sol.x;
end
sol_result(:,all(sol_result == 0,1)) = [];
sol_result_mean = mean(sol_result,2);
cd(current_path)