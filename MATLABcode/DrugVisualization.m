% Created by Taylor Eggertsen 1/31/2023

clear all;clc; close all;

model(1) = {'TAliHypertrophy_Bromocriptine_v3.xlsx'};
model(2) = {'TAliHypertrophy_Escitalopram_v4.xlsx'};
model(3) = {'TAliHypertrophy_Crizotinib_v2.xlsx'};
model(4) = {'TAliHypertrophy_PLK1combo.xlsx'};
model(5) = {'TAliHypertrophy_Tirbanibulin_v2.xlsx'};
model(6) = {'TAliHypertrophy_OSI-930_v2.xlsx'};
model(7) = {'TAliHypertrophy_Rifabutin_v2.xlsx'};
model(8) = {'TAliHypertrophy_Mifepristone_v3.xlsx'};

for mm=1:8
    warning off;
    formattedReactions = table;
    % Species/Reaction information from toy_model.xlsx or network, 'species'/'reactions' tab 
    networkReactions = readtable(model{mm}, 'Sheet', 'reactions');
    if ismember('=',cell2mat(networkReactions{1,3}))==0
        networkReactions(1,:)=[];  
    end
    % Formats network reactions to only show the product/output. 
    for i = 1:height(networkReactions)
        reaction = string(networkReactions{i,3});
        nodeOfReaction = extractAfter(reaction, '=>'); 
        formattedReactions{i,1} = strtrim(nodeOfReaction);
    end
    formattedReactions.Properties.VariableNames(1) = {'ReactionOutputNode'};

    %% Generate the drug modified ODE and parameter files from NETFLUX
    [status, result] = generateDrugODEs(model{mm});
    
    %% Create DrugsToSimulate.csv
    % % % [status, result] = generateDrugsToSimulate(model);
    
    %% Inputs for simulations
    opts = detectImportOptions('DrugsToSimulatev2.csv');
    opts = setvartype(opts,{'AgonistTargetIndex','AntagonistTargetIndex'},'char');
    drugsToSimulate = readtable('DrugsToSimulatev2.csv');% 'Z:\Taylor\Code\PPI\PPI_Edges\DrugsToSimulate2.csv',opts);% 
    
    % formattedReactions(1,:) = []; % remove empty (first) row

    %% Set drug dose or doses (as a vector) & Sensitivity analysis

    alteration_antag = [0.8] ;%%% [0:0.02:1]
    % % % [padding,list_numDrugs, list_drugAction, list_drugType, list_drugTarget, rowLabels_drugs] = createRowLabels(drugsToSimulate);
    
    inputs = ["AngII","ANPi","BNPi","CT1","EGF","ET1","FGF","IGF1","IL6","ISO","LIF","NE","NRG1","PE","Stretch","TGFB","TNFa"];
    for inp = 1:length(inputs) %input to network 
        drugofinterest = 276+mm; 
        
        % Baseline simulation
        [params,y0] = tempDrugODE_params;
        [rpar,tau,ymax,speciesNames]=params{:};
        w = rpar(1,:);
        n = rpar(2,:);
        EC50 = rpar(3,:);
        dose = rpar(4,:);
        drugBinding = rpar(5,:); 
        drugAgonism = rpar(6,:);
        w(find(ismember(networkReactions.Rule,'=>lig'))) = 0.06;
        rpar = [w;n;EC50;dose;drugBinding;drugAgonism];
        params=[rpar,tau,ymax,speciesNames];
        tspan = [0 50]; options = [];
        [tb,yb] = ode15s(@tempDrugODE,tspan,y0,options,params);
        ybase = real(yb(end,:)');
    
        % Node parameters
        [params,y0] = tempDrugODE_params;
        [rpar,tau,ymax,speciesNames]=params{:};
        w = rpar(1,:);
        n = rpar(2,:);
        EC50 = rpar(3,:);
        dose = rpar(4,:);
        drugBinding = rpar(5,:); 
        drugAgonism = rpar(6,:);
        %New input value
        w(inp) = 0.1;
        w(find(ismember(networkReactions.Rule,'=>lig'))) = 0.06;
 
        rpar = [w;n;EC50;dose;drugBinding;drugAgonism];
        inputNodeW = num2cell(1:length(speciesNames)); % Nodes to test drug against
        params=[rpar,tau,ymax,speciesNames];
        % Steady-state Control Simulation
        [t,y] = ode15s(@tempDrugODE,tspan,y0,options,params);
        yEnd = y(end,:)';
        
        % Reset initial y values
        y0 = real(yEnd); 
    
        for i = drugofinterest-1 %input drug of interest
            drugsToSimulate.Drug{i}
            drugsToSimulate.AgonistTarget{i}
            drugsToSimulate.AntagonistTarget{i}
            for j = 1:length(alteration_antag) 
                [t2,y2] = drugtreatment(dose,drugBinding,drugAgonism,drugsToSimulate,formattedReactions,alteration_antag,w,n,EC50,tau,ymax,speciesNames,y0,i,j);
                ySimEnd = y2(end,:)';
        
                specset = 19; %cellarea phenotype
                responseset(mm,inp) = (y2(end,specset)-y0(specset));
                %responseset(mm,inp) = (y2(end,specset)-y0(specset))./(y0(specset)-ybase(specset));
            end
            drugs(mm) = drugsToSimulate.Drug(i);
        end
    end
end

figure;
set(gca, 'Visible', 'on');
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
imagesc(responseset,[-1,1]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(inputs));
set(gca,'XTickLabel',inputs,'fontsize',10);
xlabel('Phenotypic Outputs','fontsize',20);
% % % xticklabel_rotate([], 45);
set(gca,'YTick',1:length(drugs));
set(gca,'YTickLabel',drugs,'fontsize',10);
ylabel('Drugs','fontsize',16);
% % % title([strcat('Change in ',{' '},specs2{dpos})]);
xtickangle(90)
hcb=colorbar;
set(get(hcb,'label'),'string','Change in Cell Area');
set(get(hcb,'label'),'fontsize',16); set(get(hcb,'label'),'rotation',90);

cgo = clustergram(responseset,'Colormap',myrgbcmap,'RowLabels',drugs,'ColumnLabels',inputs);