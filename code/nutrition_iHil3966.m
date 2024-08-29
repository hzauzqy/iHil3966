addpath ./code/nutritionAlgorithm/
load('iHil3966.mat');
parsedGPR = GPRparser_xl(model);% Extracting GPR data from the model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model); 

initCobraToolbox(false);

%disp(excludeReactions(1:10));
model2=model;

for i = 1:length(model2.rxns)
    if endsWith(model2.rxns{i}, 'e')  % Check if the reaction ends with 'e'
        model2.rxns{i} = ['EX_', model2.rxns{i}];  % Add 'EX_' at the beginning of the reaction
    end
end
% Find all reactions that start with 'EX_'
ex_reactions = strncmp('EX_', model2.rxns, 3);
% Set the upper bound of these reactions to 1
model2.ub(ex_reactions) = 1;
[model_n2] = convert_EX_to_diet(model2);

%%
solution = optimizeCbModel(model, 'max');
f1=solution.f*0.9;
model_n2.lb(ismember(model_n2.rxns,{'Biomass'})) = f1;%biomass
model_n2.lb(ismember(model_n2.rxns,{'MAR10033'})) = f1*0.08;%
model_n2.lb(ismember(model_n2.rxns,{'MAR00187'})) = f1*0.08*0.08*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model_n2.lb(ismember(model_n2.rxns,{'MAR00195'})) = f1*0.08*0.04*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model_n2.lb(ismember(model_n2.rxns,{'MAR00216'})) = f1*0.08*0.21*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model_n2.lb(ismember(model_n2.rxns,{'MAR00248'})) = f1*0.08*0.04*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model_n2.lb(ismember(model_n2.rxns,{'MAR00262'})) = f1*0.08*0.34*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model2.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.08*0.19*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model_n2.lb(ismember(model_n2.rxns,{'MAR00397'})) = f1*0.08*0.19*(279.816/279.428);
model_n2.ub(ismember(model_n2.rxns,{'MAR10033'})) = f1*0.08;%
model_n2.ub(ismember(model_n2.rxns,{'MAR00187'})) = f1*0.08*0.08*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model_n2.ub(ismember(model_n2.rxns,{'MAR00195'})) = f1*0.08*0.04*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model_n2.ub(ismember(model_n2.rxns,{'MAR00216'})) = f1*0.08*0.21*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model_n2.ub(ismember(model_n2.rxns,{'MAR00248'})) = f1*0.08*0.04*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model_n2.ub(ismember(model_n2.rxns,{'MAR00262'})) = f1*0.08*0.34*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model2.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.08*0.19*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model_n2.ub(ismember(model_n2.rxns,{'MAR00397'})) = f1*0.08*0.19*(279.816/279.428);

%%
die_rxn = {'MAR00216', 'MAR08939', 'MAR02182', 'MAR00217', 'MAR10033', 'MAR09678'};%C160
%die_rxn = {'MAR09786', 'MAR08940', 'MAR00262', 'MAR10033'};%C181
%die_rxn = {'MAR00396', 'MAR09779', 'MAR00397', 'MAR10033'};%C182
die_rxn_minmax = {'max','max','max','max','max','max'};
[newDietModel,pointsModel,roiFlux,pointsModelSln,menuChanges,detailedAnalysis] = nutritionAlgorithm2(model_n2,die_rxn,die_rxn_minmax);

%%
model_n3=newDietModel;
%
solutionBefore = optimizeCbModel(model_n3, 'max', 'one');

model_n3.ub(ismember(model_n3.rxns,{'Food_Added_EX_MAM01369e'}))=-14.916;
model_n3.lb(ismember(model_n3.rxns,{'Food_Added_EX_MAM01369e'}))=-14.916;
model_n3.ub(ismember(model_n3.rxns,{'Food_Added_EX_MAM02184e'}))=-17.469;
model_n3.lb(ismember(model_n3.rxns,{'Food_Added_EX_MAM02184e'}))=-17.469;
model_n3.ub(ismember(model_n3.rxns,{'Food_Added_EX_MAM02426e'}))=-39.389;
model_n3.lb(ismember(model_n3.rxns,{'Food_Added_EX_MAM02426e'}))=-39.389;
model_n3.ub(ismember(model_n3.rxns,{'Food_Added_EX_MAM02685e'}))=-0.3534;
model_n3.lb(ismember(model_n3.rxns,{'Food_Added_EX_MAM02685e'}))=-0.3534;

%model_n3.ub(ismember(model_n3.rxns,{'Food_Added_EX_MAM02715e'}))=-3.3466;
%model_n3.lb(ismember(model_n3.rxns,{'Food_Added_EX_MAM02715e'}))=-3.3466;
model_n3.ub(ismember(model_n3.rxns,{'Food_Added_EX_MAM02993e'}))=-30.541;
model_n3.lb(ismember(model_n3.rxns,{'Food_Added_EX_MAM02993e'}))=-30.541;

model_n3.ub(ismember(model_n3.rxns,{'Food_Added_EX_MAM03089e'}))=-12.612;
model_n3.lb(ismember(model_n3.rxns,{'Food_Added_EX_MAM03089e'}))=-12.612;
model_n3.ub(ismember(model_n3.rxns,{'Food_Added_EX_MAM03135e'}))=-65.009;
model_n3.lb(ismember(model_n3.rxns,{'Food_Added_EX_MAM03135e'}))=-65.009;


%%
%7 nutritions add before and after
model_n3.ub(ismember(model_n3.rxns,{'Diet_EX_MAM01369e'}))=0;
model_n3.lb(ismember(model_n3.rxns,{'Diet_EX_MAM01369e'}))=0;
model_n3.ub(ismember(model_n3.rxns,{'Diet_EX_MAM02426e'}))=0;%-39.389;
model_n3.lb(ismember(model_n3.rxns,{'Diet_EX_MAM02426e'}))=0;%-39.389;
model_n3.ub(ismember(model_n3.rxns,{'Diet_EX_MAM02184e'}))=0;%-17.469;
model_n3.lb(ismember(model_n3.rxns,{'Diet_EX_MAM02184e'}))=0;%-17.469;
model_n3.ub(ismember(model_n3.rxns,{'Diet_EX_MAM02685e'}))=0;%-0.3534;
%model_n3.lb(ismember(model_n3.rxns,{'Diet_EX_MAM02685e'}))=0;%-0.3534;
%model_n3.lb(ismember(model_n3.rxns,{'Diet_EX_MAM02715e'}))=0;%-3.3466;
model_n3.ub(ismember(model_n3.rxns,{'Diet_EX_MAM02993e'}))=0;%-30.541;
model_n3.lb(ismember(model_n3.rxns,{'Diet_EX_MAM02993e'}))=0;%-30.541;
model_n3.ub(ismember(model_n3.rxns,{'Diet_EX_MAM03089e'}))=0;%-12.612;
model_n3.lb(ismember(model_n3.rxns,{'Diet_EX_MAM03089e'}))=0;%-12.612;
model_n3.ub(ismember(model_n3.rxns,{'Diet_EX_MAM03135e'}))=0;%-65.009;
model_n3.lb(ismember(model_n3.rxns,{'Diet_EX_MAM03135e'}))=0;%-65.009;

% Copy the original model
tempModel2 = model_n3;
reactions = {
    'Diet_EX_MAM01369e', -14.916;
    'Diet_EX_MAM02426e', -39.389;
    'Diet_EX_MAM02184e', -17.469;
    'Diet_EX_MAM02685e', -0.3534;
    %'Diet_EX_MAM02715e', -3.3466;
    'Diet_EX_MAM02993e', -30.541;
    'Diet_EX_MAM03089e', -12.612;
    'Diet_EX_MAM03135e', -65.009
};
% Set new lower bounds for all reactions
for i = 1:length(reactions)
    tempModel2 = changeRxnBounds(tempModel2, reactions{i, 1}, reactions{i, 2}, 'l');
end

% pFBA
solutionxxx = optimizeCbModel(tempModel2, 'max', 'one');

% Set bounds of all other reactions to 0
for i = 1:length(reactions)
    tempModel2 = changeRxnBounds(tempModel2, reactions{i, 1}, 0, 'b');
end

% pFBA
solutionxxx_before = optimizeCbModel(tempModel2, 'max', 'one');

% Retrieve reaction names and flux values
rxnNames = tempModel2.rxns;
fluxesBefore = solutionxxx_before.v;
fluxesAfter = solutionxxx.v;

% Organize the data into a table
T = table(rxnNames, fluxesBefore, fluxesAfter, 'VariableNames', {'Reaction', 'FluxBefore', 'FluxAfter'});
% Save the table as an Excel file
writetable(T, './Result/nutrition7_add_beforeAndafter_flux_comparison.xlsx');


%%
%single nutrient
reactions = {
    'Diet_EX_MAM01369e', -14.916;
    'Diet_EX_MAM02426e', -39.389;
    'Diet_EX_MAM02184e', -17.469;
    'Diet_EX_MAM02685e', -0.3534;
    %'Diet_EX_MAM02715e', -3.3466;
    'Diet_EX_MAM02993e', -30.541;
    'Diet_EX_MAM03089e', -12.612;
    'Diet_EX_MAM03135e', -65.009
};

results = {};

% 
tempModel = model_n3; 
solution = optimizeCbModel(tempModel, 'max', 'one');
results{end+1, 1} = 'No change'; 
results{end, 2} = solution.f;
reactionIndex = find(strcmp(model_n3.rxns, 'MAR00031'));
results{end, 3} = solution.v(reactionIndex);

% Iterate through each reaction, modifying boundaries
for i = 1:length(reactions)
    % Create a temporary model
    tempModel = model_n3;
    
    % Set bounds for selected reactions
    tempModel = changeRxnBounds(tempModel, reactions{i, 1}, reactions{i, 2}, 'b');
    
    % Set bounds of all other reactions to 0
    for k = 1:length(reactions)
        if k ~= i
            tempModel = changeRxnBounds(tempModel, reactions{k, 1}, 0, 'b');
        end
    end
    
    % pFBA
    solution = optimizeCbModel(tempModel, 'max', 'one');
    
    % result
    results{end+1, 1} = reactions{i, 1};
    results{end, 2} = solution.f;
    reactionIndex = find(strcmp(model_n3.rxns, 'MAR00031'));
    results{end, 3} = solution.v(reactionIndex);
end
resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue', 'lipid_droplet'});
filename = './Result/pfba_results_1rxns.xlsx';
writetable(resultsTable, filename);

%%Combination of two nutrients
% Define reactions and their bounds
reactions = {
    'Diet_EX_MAM01369e', -14.916;
    'Diet_EX_MAM02426e', -39.389;
    'Diet_EX_MAM02184e', -17.469;
    'Diet_EX_MAM02685e', -0.3534;
    %'Diet_EX_MAM02715e', -3.3466;
    'Diet_EX_MAM02993e', -30.541;
    'Diet_EX_MAM03089e', -12.612;
    'Diet_EX_MAM03135e', -65.009
};
% Initialize array to store results
results = {};
% Iterate over all combinations of two reactions
for i = 1:length(reactions)-1
    for j = i+1:length(reactions)
        % Create a temporary model
        tempModel = model_n3;
        
        % Set bounds for selected reactions
        tempModel = changeRxnBounds(tempModel, reactions{i, 1}, reactions{i, 2}, 'b');
        tempModel = changeRxnBounds(tempModel, reactions{j, 1}, reactions{j, 2}, 'b');
        
        % Set bounds of all other reactions to 0
        for k = 1:length(reactions)
            if k ~= i && k ~= j
                tempModel = changeRxnBounds(tempModel, reactions{k, 1}, 0, 'b');
            end
        end
        
        % Compute pFBA
        solution = optimizeCbModel(tempModel, 'max', 'one');
        
        % Store the results
        results{end+1, 1} = {reactions{i, 1}, reactions{j, 1}};
        results{end, 2} = solution.f;
        reactionIndex2 = find(strcmp(model_n3.rxns, 'MAR00031'));
        results{end, 3} = solution.v(reactionIndex2);
    end
end
% Display the results
%resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue'});
resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue', 'lipid_droplet'});
filename = './Result/pfba_results_2rxns.xlsx';
writetable(resultsTable, filename);

 %%Combination of three nutrients
 % Iterate over all combinations of three reactions
 results = {};
for i = 1:length(reactions)-2
    for j = i+1:length(reactions)-1
        for k = j+1:length(reactions)
            % Create a temporary model
            tempModel = model_n3;
            
            % Set bounds for selected reactions
            tempModel = changeRxnBounds(tempModel, reactions{i, 1}, reactions{i, 2}, 'b');
            tempModel = changeRxnBounds(tempModel, reactions{j, 1}, reactions{j, 2}, 'b');
            tempModel = changeRxnBounds(tempModel, reactions{k, 1}, reactions{k, 2}, 'b');
            
            % Set bounds of all other reactions to 0
            for m = 1:length(reactions)
                if m ~= i && m ~= j && m ~= k
                    tempModel = changeRxnBounds(tempModel, reactions{m, 1}, 0, 'b');
                end
            end
            
            % Compute pFBA
            solution = optimizeCbModel(tempModel, 'max', 'one');
            
            % Store the results
            results{end+1, 1} = {reactions{i, 1}, reactions{j, 1}, reactions{k, 1}};
            results{end, 2} = solution.f;
            reactionIndex2 = find(strcmp(model_n3.rxns, 'MAR00031'));
            results{end, 3} = solution.v(reactionIndex2);
        end
    end
end
%resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue'});
resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue', 'lipid_droplet'});
filename = './Result/pfba_results_3rxns.xlsx';
writetable(resultsTable, filename);

%%Combination of four nutrients
 results = {};
% Iterate over all combinations of four reactions
for i = 1:length(reactions)-3
    for j = i+1:length(reactions)-2
        for k = j+1:length(reactions)-1
            for l = k+1:length(reactions)
                % Create a temporary model
                tempModel = model_n3;
                
                % Set bounds for selected reactions
                tempModel = changeRxnBounds(tempModel, reactions{i, 1}, reactions{i, 2}, 'b');
                tempModel = changeRxnBounds(tempModel, reactions{j, 1}, reactions{j, 2}, 'b');
                tempModel = changeRxnBounds(tempModel, reactions{k, 1}, reactions{k, 2}, 'b');
                tempModel = changeRxnBounds(tempModel, reactions{l, 1}, reactions{l, 2}, 'b');
                
                % Set bounds of all other reactions to 0
                for m = 1:length(reactions)
                    if m ~= i && m ~= j && m ~= k && m ~= l
                        tempModel = changeRxnBounds(tempModel, reactions{m, 1}, 0, 'b');
                    end
                end
                
                % Compute pFBA
                solution = optimizeCbModel(tempModel, 'max', 'one');
                
                % Store the results
                results{end+1, 1} = {reactions{i, 1}, reactions{j, 1}, reactions{k, 1}, reactions{l, 1}};
                results{end, 2} = solution.f;
                reactionIndex2 = find(strcmp(model_n3.rxns, 'MAR00031'));
                results{end, 3} = solution.v(reactionIndex2);
            end
        end
    end
end
%resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue'});
resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue', 'lipid_droplet'});
filename = './Result/pfba_results_4rxns.xlsx';
writetable(resultsTable, filename);

%%Combination of five nutrients
 results = {};
% Iterate over all combinations of five reactions
for i = 1:length(reactions)-4
    for j = i+1:length(reactions)-3
        for k = j+1:length(reactions)-2
            for l = k+1:length(reactions)-1
                for m = l+1:length(reactions)
                    % Create a temporary model
                    tempModel = model_n3;
                    
                    % Set bounds for selected reactions
                    tempModel = changeRxnBounds(tempModel, reactions{i, 1}, reactions{i, 2}, 'b');
                    tempModel = changeRxnBounds(tempModel, reactions{j, 1}, reactions{j, 2}, 'b');
                    tempModel = changeRxnBounds(tempModel, reactions{k, 1}, reactions{k, 2}, 'b');
                    tempModel = changeRxnBounds(tempModel, reactions{l, 1}, reactions{l, 2}, 'b');
                    tempModel = changeRxnBounds(tempModel, reactions{m, 1}, reactions{m, 2}, 'b');
                    
                    % Set bounds of all other reactions to 0
                    for n = 1:length(reactions)
                        if n ~= i && n ~= j && n ~= k && n ~= l && n ~= m
                            tempModel = changeRxnBounds(tempModel, reactions{n, 1}, 0, 'b');
                        end
                    end
                    
                    % Compute pFBA
                    solution = optimizeCbModel(tempModel, 'max', 'one');
                    
                    % Store the results
                    results{end+1, 1} = {reactions{i, 1}, reactions{j, 1}, reactions{k, 1}, reactions{l, 1}, reactions{m, 1}};
                    results{end, 2} = solution.f;
                    reactionIndex2 = find(strcmp(model_n3.rxns, 'MAR00031'));
                    results{end, 3} = solution.v(reactionIndex2);
                end
            end
        end
    end
end
%resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue'});
resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue', 'lipid_droplet'});
filename = './Result/pfba_results_5rxns.xlsx';
writetable(resultsTable, filename);

%%Combination of six nutrients
 results = {};
% Iterate over all combinations of six reactions
for i = 1:length(reactions)-5
    for j = i+1:length(reactions)-4
        for k = j+1:length(reactions)-3
            for l = k+1:length(reactions)-2
                for m = l+1:length(reactions)-1
                    for n = m+1:length(reactions)
                        % Create a temporary model
                        tempModel = model_n3;
                        
                        % Set bounds for selected reactions
                        tempModel = changeRxnBounds(tempModel, reactions{i, 1}, reactions{i, 2}, 'b');
                        tempModel = changeRxnBounds(tempModel, reactions{j, 1}, reactions{j, 2}, 'b');
                        tempModel = changeRxnBounds(tempModel, reactions{k, 1}, reactions{k, 2}, 'b');
                        tempModel = changeRxnBounds(tempModel, reactions{l, 1}, reactions{l, 2}, 'b');
                        tempModel = changeRxnBounds(tempModel, reactions{m, 1}, reactions{m, 2}, 'b');
                        tempModel = changeRxnBounds(tempModel, reactions{n, 1}, reactions{n, 2}, 'b');
                        
                        % Set bounds of all other reactions to 0
                        for o = 1:length(reactions)
                            if o ~= i && o ~= j && o ~= k && o ~= l && o ~= m && o ~= n
                                tempModel = changeRxnBounds(tempModel, reactions{o, 1}, 0, 'b');
                            end
                        end
                        
                        % Compute pFBA
                        solution = optimizeCbModel(tempModel, 'max', 'one');
                        
                        % Store the results
                        results{end+1, 1} = {reactions{i, 1}, reactions{j, 1}, reactions{k, 1}, reactions{l, 1}, reactions{m, 1}, reactions{n, 1}};
                        results{end, 2} = solution.f;
                        reactionIndex2 = find(strcmp(model_n3.rxns, 'MAR00031'));
                        results{end, 3} = solution.v(reactionIndex2);
                    end
                end
            end
        end
    end
end
%resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue'});
resultsTable = cell2table(results, 'VariableNames', {'ReactionPair', 'ObjectiveValue', 'lipid_droplet'});
filename = './Result/pfba_results_6rxns.xlsx';
writetable(resultsTable, filename);
