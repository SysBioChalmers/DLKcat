function [model,sol_result_mean] = getflux(strain,chemostat,media,anox,carbon,inputpath,state)

if ~strcmp(chemostat,'max')
    dilutionrate = str2double(chemostat);
    chemostat = 'chemo';
elseif strcmp(chemostat,'max')
    chemostat = 'batch';
end
current_path = pwd;
cd(inputpath)

if strcmp(state,'prior') || strcmp(state,'auto') || strcmp(state,'median')
    num = 1;
elseif strcmp(state,'posterior')
    num = 5;
end
for m = 1:num
    m
    if strcmp(state,'prior')
        load(['ecmodel_',strain,'_prior.mat'])
    elseif strcmp(state,'posterior')
        load(['ecmodel_',strain,num2str(m),'.mat'])
    elseif strcmp(state,'median')
        cd('../model_auto')
        z = load([strain,'.mat']);
        if strcmp(strain,'Kluyveromyces_marxianus')
            tot_prot = 0.325;
        elseif strcmp(strain,'Kluyveromyces_lactis')
             tot_prot = 0.245;
        else
             tot_prot = 0.23;
        end
        z.enzymedata.kcat(1:length(z.enzymedata.kcat)) = median(z.enzymedata.kcat);
        emodel = convertToGeckoModel(z.model,z.enzymedata,tot_prot);
    else
        cd('../model_auto')
        z = load([strain,'.mat']);
        if strcmp(strain,'Kluyveromyces_marxianus')
            tot_prot = 0.325;
        elseif strcmp(strain,'Kluyveromyces_lactis')
             tot_prot = 0.245;
        else
             tot_prot = 0.23;
        end
        z.enzymedata.kcat(z.enzymedata.kcat < 36) = 36;
        z.enzymedata.kcat(z.enzymedata.kcat > 36000000000) = 36000000000;
        emodel = convertToGeckoModel(z.model,z.enzymedata,tot_prot);
        cd ../
        cd(inputpath)
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