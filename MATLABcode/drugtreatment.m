function [t2,y2,targetnode] = drugtreatment(dose,drugBinding,drugAgonism,drugsToSimulate,formattedReactions,alteration_antag,w,n,EC50,tau,ymax,speciesNames,y0,i,j)
% Written by Anirudha Chandrabhatla 
% Modified by Taylor Eggertsen 
% Last update 7/12/2022
% Simulates the acivity of the drug in the network model according to its
% dose, as well as binding and agonistic properties.

doseNew = dose;
drugBindingNew = drugBinding;
drugAgonismNew = drugAgonism;

% The drug is an agonist
if strcmp(drugsToSimulate.IsAgonist{i}, 'Yes') == 1
    if isempty(find(drugsToSimulate.AgonistTarget{i} == ';', 1)) % Drug has one agonist target
        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTarget{i});
        targetnode = drugsToSimulate.AgonistTarget{i};
        if strcmp('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, agonist
            drugBindingNew(locationOfReactions) = 1;
            drugAgonismNew(locationOfReactions) = 1;
            doseNew(locationOfReactions) = -1*alteration_antag(j);
        else % Non-Competitive, agonist
            drugBindingNew(locationOfReactions) = -1;
            drugAgonismNew(locationOfReactions) = 1;
            doseNew(locationOfReactions) = alteration_antag(j);
        end
    else % Drug has multiple targets
        geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTarget{i}, ';');
        targetnode = geneIDsOfTargets;
        for m = 1:length(geneIDsOfTargets)
            locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{m});
            if strcmp('Competitive', drugsToSimulate.DrugAction{i})
                drugBindingNew(locationOfReactions) = 1;
                drugAgonismNew(locationOfReactions) = 1;
                doseNew(locationOfReactions) = -1*alteration_antag(j);
            else 
                drugBindingNew(locationOfReactions) = -1;
                drugAgonismNew(locationOfReactions) = 1;
                doseNew(locationOfReactions) = alteration_antag(j);
            end
        end
    end
end

% The drug is an antagonist
if strcmp(drugsToSimulate.IsAntagonist{i}, 'Yes') == 1
    if isempty(find(drugsToSimulate.AntagonistTarget{i} == ';', 1)) % Drug has one antagonist target
        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AntagonistTarget{i});
        targetnode = drugsToSimulate.AntagonistTarget{i};
        if strcmp('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, antagonist
            drugBindingNew(locationOfReactions) = 1;
            drugAgonismNew(locationOfReactions) = -1;
            doseNew(locationOfReactions) = alteration_antag(j);
        else
            drugBindingNew(locationOfReactions) = -1;
            drugAgonismNew(locationOfReactions) = -1;
            doseNew(locationOfReactions) = alteration_antag(j);
        end
    else
        geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTarget{i}, ';');
        targetnode = geneIDsOfTargets;
        for p = 1:length(geneIDsOfTargets)
            locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
            if strcmp('Competitive', drugsToSimulate.DrugAction{i})
                drugBindingNew(locationOfReactions) = 1;
                drugAgonismNew(locationOfReactions) = -1;
                doseNew(locationOfReactions) = alteration_antag(j);
            else 
                drugBindingNew(locationOfReactions) = -1;
                drugAgonismNew(locationOfReactions) = -1;
                doseNew(locationOfReactions) = alteration_antag(j);
            end
        end
    end
end

rparNew = [w;n;EC50;doseNew;drugBindingNew;drugAgonismNew];
paramsNew = {rparNew,tau,ymax,speciesNames};
tspan = [0 50]; options = []; 
[t2,y2] = ode15s(@tempDrugODE,tspan,y0,options,paramsNew); % Make sure y0 is correct here.
end

