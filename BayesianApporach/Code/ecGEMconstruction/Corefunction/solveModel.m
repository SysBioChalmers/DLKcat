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
    sol.exitflag = resultgurobi.status;
    if strcmp(resultgurobi.status,'OPTIMAL')
        sol.x = resultgurobi.x;
        sol.obj = sol.x(logical(model.c));
        sol.protUsage = cost_list * sol.x / 1000;
    else
        sol.obj = [];
        sol.protUsage = [];
        sol.x = zeros(length(model.rxns),1);
    end
    
elseif  strcmp(solver,'ibm_cplex')
            param.verify= 0;
       param.minNorm= [];
    param.printLevel= 0;
    param.primalOnly= 0;
     param.saveInput= [];
       param.feasTol= 1.0000e-06;
        param.optTol= 1.0000e-06;
        param.solver= 'ibm_cplex';
         param.debug= 0;
       param.logFile= [];
       param.lifting= 0;
        param.method= -1;
        solverParams = struct();
    % set constraints
    optProblem.A = [model.S;cost_list];
    optProblem.lb = model.lb;
    optProblem.ub = model.ub;
    optProblem.b = [zeros(nMets,1);tot_prot_weight*1000];
    optProblem.csense = [repmat('E',nMets,1);'L'];
    if strcmp(osenseStr,'max')
        optProblem.osense = -1;
    else
        optProblem.osense = 1;
    end
    optProblem.c = model.c;
    % call the solver
    CplexProblem = buildCplexProblemFromCOBRAStruct(optProblem); % a cobra function to change from the constraint to cplex format
            [CplexLPproblem, logFile, logToFile] = setCplexParametersForProblem(CplexProblem,param,solverParams,'LP');
        %logToFile=0;
        
        % optimize the problem
        CplexLPproblem.solve();
    
        if logToFile
            % Close the output file
            fclose(logFile);
        end
                
        % http://www-eio.upc.edu/lceio/manuals/cplex-11/html/overviewcplex/statuscodes.html
        % https://www.ibm.com/support/knowledgecenter/SSSA5P_12.5.1/ilog.odms.cplex.help/refmatlabcplex/html/classCplex.html#a93e3891009533aaefce016703acb30d4
        origStat   = CplexLPproblem.Solution.status;
        %stat = origStat;
        if origStat==1 && isfield(CplexLPproblem.Solution,'dual')
            stat = 1;
            f = CplexLPproblem.Solution.objval;
            if ~isfield(CplexLPproblem.Solution,'x')
                disp(CplexLPproblem)
            end
            sol.x = CplexLPproblem.Solution.x;
            sol.obj = sol.x(logical(model.c));
            sol.protUsage = cost_list * sol.x / 1000;
        elseif origStat == 2 ||   origStat == 20
            stat = 2; %unbounded
             sol.obj = [];
        sol.protUsage = [];
        sol.x = zeros(length(model.rxns),1);
        
        elseif origStat == 3
            stat = 0;%infeasible
             sol.obj = [];
        sol.protUsage = [];
        sol.x = zeros(length(model.rxns),1);
        elseif origStat == 4
            % this is likely unbounded, but could be infeasible
            % lets check, by solving an additional LP with no objective.
            % if that LP has a solution, it's unbounded. If it doesn't, it's infeasible.
            Solution = CplexLPproblem.Solution;
            CplexLPproblem.Model.obj(:) = 0;
            CplexLPproblem.solve();
            origStatNew   = CplexLPproblem.Solution.status;
            if origStatNew == 1
                stat = 2;
                 sol.obj = [];
        sol.protUsage = [];
        sol.x = zeros(length(model.rxns),1);
            else
                stat = 0;
                 sol.obj = [];
        sol.protUsage = [];
        sol.x = zeros(length(model.rxns),1);
            end
            % restore the original solution.
            % restore the original solution.
            CplexLPproblem.Solution = Solution;
        elseif origStat == 5 || origStat == 6
            stat = 3;% Almost optimal solution
            sol.x = CplexLPproblem.Solution.x;
            sol.obj = sol.x(logical(model.c));
            sol.protUsage = cost_list * sol.x / 1000;
        elseif (origStat >= 10 && origStat <= 12) || origStat == 21 || origStat == 22
            % abort due to reached limit. check if there is a solution and return it.
            stat = 3;
            if isfield(CplexLPproblem.Solution ,'x')
                sol.x = CplexLPproblem.Solution.x;
            else
               % no solution returned
                stat = -1;
            end
        else
            stat = -1;
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