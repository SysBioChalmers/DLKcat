function sol = solveModel(model,objective,osenseStr,prot_cost_info,tot_prot_weight,solver,rxnlst,factor)



% Determine objective
model.c(strcmp(model.rxns,objective)) = 1;
changerxn = true;
if nargin < 7
    changerxn = false;
end
% Add protein cost infomation
[nMets,nRxns] = size(model.S);
cost_list = zeros(1,nRxns);
for i = 1:nRxns
    rxnid = model.rxns{i};
    if changerxn && any(strcmp(rxnlst,rxnid))
        id_tmp = rxnid;
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp))/factor;
        else
            cost = 0;
        end
    else
        id_tmp = rxnid;
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp));
        else
            cost = 0;
        end
    end
    cost_list(1,i) = cost;
end
if strcmp(solver,'gurobi')
    param.DisplayInterval = 1;
    param.FeasibilityTol = 1.0000e-06;
    param.OptimalityTol = 1.0000e-06;
    param.OutputFlag = 0;
    % set constraints
    gurobiLP.A = [model.S;cost_list];
    gurobiLP.lb = model.lb;
    gurobiLP.ub = model.ub;
    gurobiLP.rhs = [zeros(nMets,1);tot_prot_weight*1000];
    gurobiLP.sense = [repmat('=',nMets,1);'<'];
    gurobiLP.modelsense = osenseStr;
    gurobiLP.obj = model.c;
    % call the solver
    resultgurobi = gurobi(gurobiLP,param);
    sol = struct();
    sol.x = resultgurobi.x;
    sol.exitflag = resultgurobi.status;
    if strcmp(resultgurobi.status,'OPTIMAL')
        sol.obj = sol.x(logical(model.c));
        sol.protUsage = cost_list * sol.x / 1000;
    else
        sol.obj = [];
        sol.protUsage = [];
    end
else % use matlab solver
    % Construct LP
    A = cost_list;
    b = tot_prot_weight*1000;
    
    Aeq = model.S;
    beq = zeros(nMets,1);
    
    lb = model.lb;
    ub = model.ub;
    
    %linprog always runs minimization
    if osenseStr == 'max'
        f = -model.c;
    elseif osenseStr == 'min'
        f = model.c;
    end
    
    [x,~,exitflag,~] = linprog(f,A,b,Aeq,beq,lb,ub);
    
    sol = struct();
    sol.x = x;
    sol.exitflag = exitflag;
    if exitflag == 1
        sol.obj = sol.x(logical(model.c));
        sol.protUsage = cost_list * sol.x / 1000;
    else
        sol.obj = [];
        sol.protUsage = [];
    end

end