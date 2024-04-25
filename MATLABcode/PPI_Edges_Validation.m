%% Identify connections using PathLinker
% created by Taylor Eggertsen
% last updated 11/15/2022

clear all; clc; close all;
AddGenes = readtable('DrugGenes.csv','ReadVariableNames',true);
model='Hypertrophy-Model.xlsx';

%Read in network model and pull network genes out
network = readtable(model,'sheet','species');
nodenames = table2array(network(:,2));
genenames = table2array(network(:,8));

% import data from PathLinker
allnodes = readtable('SignalingNetworks_Human_KRS.csv');
edges = readtable('PathLinker-network.csv'); %change file source
node1 = strtrim(string(edges.source_genesymbol));
node2 = strtrim(string(edges.target_genesymbol));
rank = edges.pathRank4; %change appropriate column  name
con = edges.is_inhibition;
chk = edges.is_stimulation;
for i = 1:length(con)
    if (con(i)+chk(i)) == 2
        con(i) = 2;
        
    end
end

% Create new table displaying reactions, listed by rank
rel = [];
for i = 1:length(rank)
    if con(i)==1 && chk(i)==0
        rel = vertcat(rel,"Inhibits");
    elseif con(i)==0 && chk(i)==1
        rel = vertcat(rel,"Stimulates");
    else
        rel = vertcat(rel,"Stimulates/Inhibits");
    end
end
Connections = table(node1,rel,node2,rank);
Connections_rank = table2array(sortrows(Connections,4));
Connections_name = table2array(sortrows(Connections,1));

%Rewrite edges in the Netflux format to be easily integratable
A=[]; An=[]; B=[]; Bn=[]; C=[]; 
for i = 1:length(Connections_rank)
    if Connections_rank(i,2)=="Stimulates"
        A = vertcat(A,Connections_rank(i,1));
        B = vertcat(B,Connections_rank(i,3));
        C = vertcat(C,Connections_rank(i,4));
    elseif Connections_rank(i,2)=="Inhibits"
        A = vertcat(A,strcat("!",Connections_rank(i,1)));
        B = vertcat(B,Connections_rank(i,3));
        C = vertcat(C,Connections_rank(i,4));
    elseif Connections_rank(i,2)=="Stimulates/Inhibits"
        A = vertcat(A,Connections_rank(i,1));
        B = vertcat(B,Connections_rank(i,3));
        C = vertcat(C,Connections_rank(i,4));
        A = vertcat(A,strcat("!",Connections_rank(i,1)));
        B = vertcat(B,Connections_rank(i,3));
        C = vertcat(C,-1*str2double(Connections_rank(i,4)));
    end
end
for i = 1:length(B)
    index = contains(genenames,B(i));
    if sum(index)==1
        Bn = vertcat(Bn,nodenames(find(index)));
    elseif sum(index)>1
        pos = find(index); 
        for j=1:length(pos)
            set = strsplit(genenames{pos(j)},';');
            count = 0;
            for k=1:length(set)
                count = count + sum(strmatch(set{k},B(i)));
            end
            if count>0
                Bn = vertcat(Bn,nodenames(pos(j)));
            end
        end
    elseif sum(index)==0
        Bn = vertcat(Bn,B(i));
    end
end
for i = 1:length(A)
    index = contains(genenames,erase(A(i),'!'));
    if sum(index)==1
        if contains(A(i),'!')==1
            An = vertcat(An,strcat('!',nodenames(find(index))));
        else
            An = vertcat(An,nodenames(find(index)));
        end
    elseif sum(index)>1
        pos = find(index); 
        for j=1:length(pos)
            set = strsplit(genenames{pos(j)},';');
            count = 0;
            for k=1:length(set)
                count = count + sum(strmatch(set{k},A(i)));
            end
            if count>0
                An = vertcat(An,nodenames(pos(j)));
            end
        end
    elseif sum(index)==0
        An = vertcat(An,A(i));
    end
end
EdgesFormatted = [strcat(An," => ",Bn),C];
[~,uid]=unique(EdgesFormatted(:,1),'stable');
An = An(uid); Bn = Bn(uid); C = C(uid);
EdgesFormatted = [strcat(An," => ",Bn),C];

% multi-step pathways represented
interm=[];
for i = 1:length(An)
    if isempty(find(contains(AddGenes.Targets,erase(An(i),'!')), 1))==1 && isempty(find(contains(nodenames,erase(An(i),'!')), 1))==1
        interm = vertcat(interm,erase(An(i),'!'));
    end
end
interm = unique(interm);
for i=1:length(interm)
    multindex = find(contains(An,interm(i)));
    if length(multindex)>1
        for k=1:length(multindex)
            if isempty(multindex(k))==0 && str2num(C(multindex(k)))>0
                ranks = C(find(contains(Bn,interm(i))));
                for j=1:length(ranks)
                    EdgesFormatted = vertcat(EdgesFormatted,[strcat(An(multindex(k))," => ",Bn(multindex(k))),ranks(j)]);
                end
            end
        end
    else
        if isempty(multindex)==0
            ranks = C(find(contains(Bn,interm(i))));
            for j=1:length(ranks)
                EdgesFormatted = vertcat(EdgesFormatted,[strcat(An(multindex)," => ",Bn(multindex)),ranks(j)]);
            end
        end
    end
end

perc = ["Gene","Percent Represented in Paths"];
for i=1:height(AddGenes)
    rep = strfind(EdgesFormatted(:,1),AddGenes.Targets{i});
    lrep = length(find(not(cellfun('isempty',rep))));
    percnum = 100*lrep/length(EdgesFormatted);
    perc = vertcat(perc,[AddGenes.Targets{i},string(percnum)]);
end

nothere = ["Genes NOT in Omnipath"];
for i=1:height(AddGenes)
    rep = strfind(allnodes.gene_symbol,AddGenes.Targets{i});
    lrep = length(find(not(cellfun('isempty',rep))));
    if lrep==0
        nothere = vertcat(nothere,AddGenes.Targets{i});
    end
end

[percset,ploc] = setdiff(perc(2:end,1),nothere,'stable','rows');

p = categorical(percset);
q = str2double(perc(2:end,2));
q = q(ploc);
figure
bar(p,q)
xlabel("New Genes")
ylabel("Genes in Identified Pathways (%)")

%% Test new pathways against validation sheet
GSValidationfile_path ='Gold_Standard_Validation.xlsx';
EdgeValidation = [];
Int_time = 40;
Steady_time = 40;
Threshold = 0.1; % percent of changes
Model_version = 2; % 1= Original 2= Modified
[GpercentMatch, GresultChart, GBMatch, GbyClass] = Automated_Validation_V1(model, GSValidationfile_path, Int_time, Steady_time,Threshold, Model_version);
EdgeValidation = vertcat(EdgeValidation,[0,GpercentMatch,0,0,0,0,0]);
stiminhib = 1; %choose stimulate (1), choose inhibit (2)
winp = 0.1; %weight of input reactions

for i=min(str2double(EdgesFormatted(:,2))):max(str2double(EdgesFormatted(:,2)))
    disp(i)
    % select the the edge (or paired edges) to add to the file
    index = find(str2double(EdgesFormatted(:,2))==i);
    if isempty(index)==0
        newedge = EdgesFormatted(index,1);
        newedge = unique(newedge);
    
        inputFolder = pwd;
        sourceFile = fullfile(inputFolder, 'TAliHypertrophy_0.xlsx');  
        outputFolder = inputFolder;
        %delete previous pathlinker files
        d = dir('*.xlsx'); 
        for j = 1:length(d)
            if contains(d(j).name,'_PL_')
                delete(d(j).name)  
            end
        end
        
        outputFile = fullfile(outputFolder,'TAliHypertrophy_PL.xlsx');
        copyfile(sourceFile, outputFile);
        networkOutputs = table; edgeOutputs = table; edgeInputs = table;
        % Species/Reaction information from network, 'species'/'reactions' tab 
        networkReactions = readtable(model, 'Sheet', 'reactions');
        if ismember('=',cell2mat(networkReactions{1,3}))==0
            networkReactions(1,:)=[];  
        end
        % Formats network reactions to only show the product/output. 
        for j = 1:height(networkReactions)
            reaction = string(networkReactions{j,3});
            nodeOfReaction = extractAfter(reaction, '=>');
            networkOutputs{j,1} = strtrim(nodeOfReaction);
        end
        addedges = '';
        for j = 1:height(newedge)
            edgeOutputs{j,1} = strtrim(extractAfter(string(newedge(j)),'=>'));
            edgeInputs{j,1} = strtrim(extractBefore(string(newedge(j)),'=>'));
            % identify the new input
            count1 = strfind(table2array(edgeOutputs),erase(edgeInputs{j,1},'!'));
            count2 = strfind(table2array(networkOutputs),erase(edgeInputs{j,1},'!'));
            if iscell(count1)==0
                count1 = num2cell(count1);
            end
            if sum(cell2mat(count1))+sum(cell2mat(count2))==0
                newinput = strcat('''=>',erase(edgeInputs{j,1},'!'));
            end
            newspecies(j) = erase(edgeInputs{j,1},'!');
            addedges = addedges + newedge(j);
            if j<height(newedge)
                addedges = addedges + ';';
            end
        end
        if isempty(find(contains(nodenames,edgeInputs.Var1), 1))==1
            
            % create a new model file including the new edge(s)
            for j = 1:length(newedge)
                [~, txt,~] = xlsread('TAliHypertrophy_PL.xlsx','reactions');
                if strfind(newedge(j),'!')>0 & isempty(index) == 0
                    loc = find(table2array(networkOutputs) == table2array(edgeOutputs(j,:)));
                    for k=1:length(loc)
                        merged = strcat(table2array(edgeInputs(j,:))," & ",table2array(networkReactions(loc(k),3)));
                        xlswrite('TAliHypertrophy_PL.xlsx', ...
                            [{'middle'},{strcat('r',num2str(loc(k)))},merged,1,1.223,0.5],'reactions',strcat('A',num2str(loc(k)+2)));
                    end
                else
                    xlswrite('TAliHypertrophy_PL.xlsx', ...
                                    [{'middle'},{strcat('r',num2str(size(txt,1)+1))},newedge(j),1,1.223,0.5],'reactions',strcat('A',num2str(size(txt,1)+1)));
                end
            end
            [~, txt,~] = xlsread('TAliHypertrophy_PL.xlsx','reactions');
            [~, txt2,~] = xlsread('TAliHypertrophy_PL.xlsx','species');
            xlswrite('TAliHypertrophy_PL.xlsx', ...
                                    [{'middle'},{strcat('r',num2str(size(txt,1)+1))},newinput,winp,1.223,0.5],'reactions',strcat('A',num2str(size(txt,1)+1)));
            for k=1:length(newspecies)
                xlswrite('TAliHypertrophy_PL.xlsx', ...
                                    [{'macromolecule'},newspecies(k),newspecies(k),0,1,1,'protein'],'species',strcat('A',num2str(size(txt2,1)+k)));
            end
        
            Modelfile_path ='TAliHypertrophy_PL.xlsx';
        
            % Test the effect of the new input on cell area
                % Generate the ODE and parameter files from NETFLUX
            if exist([pwd '\ODEfun.m'],'file') == 2
                delete('ODEfun.m');
            end
            if exist([pwd '\ODEfun_loadParams.m'],'file') == 2
                delete('ODEfun_loadParams.m');
            end
            namepos = findstr('.xls', Modelfile_path); namestr = cellstr(Modelfile_path(1:namepos-1));
            
            [specID,reactionIDs,~,paramList,ODElist,~, error] = util.xls2Netflux(namestr,Modelfile_path);
            commandODE = util.exportODE2(specID,paramList,ODElist);
            [a,commandPARAM,b] = util.exportODE(specID,paramList,ODElist,'ODEfun');
            util.textwrite('ODEfun.m',commandODE);
            util.textwrite('ODEfun_loadParams.m',commandPARAM);
        
            [params,y0] = ODEfun_loadParams;
            [rpar,tau,ymax,speciesNames]=params{:};
            w = rpar(1,:);
            w(14) = 0.1;
        
            n = rpar(2,:);
            EC50 = rpar(3,:);
            rpar = [w;n;EC50];
            params=[rpar,tau,ymax,speciesNames];
            % Steady-state Control Simulation
            tspan = [0 50]; options = [];
            [t,y] = ode15s(@ODEfun,tspan,y0,options,params);
            % Reset initial y values
            y0 = real(y(end,:)');
        
            ymax(find(contains(speciesNames,strtrim(extractAfter(string(newinput),'=>'))))) = 0;
            params=[rpar,tau,ymax,speciesNames];
            tspan = [0 50]; options = [];
            [t2,y2] = ode15s(@ODEfun,tspan,y0,options,params);
            ySimEnd = real(y2(end,:)');
            activityChange = ySimEnd - y0;
        
            % run that against the validation!
            Int_time = 40;
            Steady_time = 40;
            Threshold = 0.1; % percent of changes
            Model_version = 2; % 1= Original 2= Modified
            if abs(activityChange(19))>0.001
                [GpercentMatch, GresultChart, GBMatch, GbyClass] = Automated_Validation_V1(Modelfile_path, GSValidationfile_path, Int_time, Steady_time,Threshold, Model_version);
                GpercentMatch   
            else
                GpercentMatch = "N/A";
            end
        
            if isempty(index)==0
                EdgeValidation = vertcat(EdgeValidation,[addedges,GpercentMatch,activityChange(19),activityChange(12),activityChange(4),activityChange(75),activityChange(72)]);
            end
        end
    end
end
clear set;
[z,iz,ic] = unique(EdgeValidation(:,1));
EdgeValidation = EdgeValidation(iz,:);
save('EdgeValidation.mat','EdgeValidation');

% add drug name column
drugs = cell(length(EdgeValidation),1);
label = cell(length(EdgeValidation),1);
for i = 1:height(AddGenes)
    drugx = contains(EdgeValidation(:,1),AddGenes.Targets{i});
    drugs(drugx) = AddGenes.Drug(i);
end
for i = 1:height(EdgeValidation)
    label(i) = {strcat(drugs{i}," | ",EdgeValidation{i,1})};
end

% % % % remove poorly validating edges
% % % thresh = 78; newindex = [];
% % % for i = 1:length(EdgeValidation)
% % %     if str2double(EdgeValidation(i,2)) > thresh && abs(str2double(EdgeValidation(i,3)))>0.001
% % %         newindex = vertcat(newindex,i);
% % %     end
% % % end
% % % EdgeValidation = EdgeValidation(newindex,:);
% % % save('EdgeValidation_small.mat','EdgeValidation');

% remove pathways that have no impact on hypertrophy
newindex = [];
for i = 1:length(EdgeValidation)
    if abs(str2double(EdgeValidation(i,3))) > (10^-2)
        newindex = vertcat(newindex,i);
    end
end
EdgeValidation = EdgeValidation(newindex,:);
label = label(newindex);
EdgeValidation(:,1) = label;

evnames = EdgeValidation(:,1);
x = categorical(evnames);
% % % x = reordercats(x,EdgeValidation(:,1));
yarea = str2double(EdgeValidation(:,3));
z = str2double(EdgeValidation(:,2));
ybnp = str2double(EdgeValidation(:,4));
yamhc = str2double(EdgeValidation(:,5));
ynfkb = str2double(EdgeValidation(:,6));
ymtor = str2double(EdgeValidation(:,7));
[sortedY, sortOrder] = sort(EdgeValidation(:,3),'ascend');
sortedx = x(sortOrder); 
sortedx = reordercats(sortedx,evnames(sortOrder));
sortedz = z(sortOrder);
sortedybnp = ybnp(sortOrder);
sortedyamhc = yamhc(sortOrder);
sortedynfkb = ynfkb(sortOrder);
sortedymtor = ymtor(sortOrder);
figure
barh(sortedx,str2double(sortedY));
xlabel('Change in Cell Area')

% % % figure
% % % hAx1 = subplot(1,3,1);
% % % pos = get(hAx1,'Position'); pos(1)=0.2; pos(3)=0.2;
% % % set(hAx1,'Position',pos)
% % % barh(sortedx,str2double(sortedY));
% % % xlabel('Change in Cell Area')
% % % 
% % % hAx2 = subplot(1,3,2);
% % % pos = get(hAx2,'Position'); pos(1)=0.425; pos(3)=0.2;
% % % set(hAx2,'Position',pos)
% % % barh(sortedx,sortedzz)
% % % set(gca,'YTickLabel',{' '})
% % % xlabel('Change in ANP')
% % % 
% % % hAx3 = subplot(1,3,3);
% % % pos = get(hAx3,'Position'); pos(1)=0.65; pos(3)=0.2;
% % % set(hAx3,'Position',pos)
% % % barh(sortedx,sortedz)
% % % xline(80,'r--')
% % % set(gca,'YTickLabel',{' '})
% % % xlim([75 85])
% % % xlabel('Validation Percentage')

heat = str2double(horzcat(sortedY,sortedybnp,sortedyamhc,sortedynfkb,sortedymtor)); 
heat = flipud(heat); clear set;

figure
hAx1 = subplot(1,2,1);
pos = get(hAx1,'Position'); pos(1)=0.3; pos(3)=0.35;
set(hAx1,'Position',pos)
set(gca, 'Visible', 'on');
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
imagesc(heat,[-1,1]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:5);
set(gca,'XTickLabel',["Cell Area","BNP","aMHC","NFkB","MTOR"],'fontsize',10);
xlabel('Phenotypic Outputs','fontsize',20);
set(gca,'YTick',1:length(sortedx));
set(gca,'YTickLabel',sortedx,'fontsize',10);
% ylabel('Representative Drugs','fontsize',16);
xtickangle(0)
hcb=colorbar;
set(get(hcb,'label'),'string','Change in Activity');
set(get(hcb,'label'),'fontsize',16); set(get(hcb,'label'),'rotation',90);

hAx3 = subplot(1,2,2);
pos = get(hAx3,'Position'); pos(1)=0.75; pos(3)=0.15;
set(hAx3,'Position',pos)
barh(sortedx,sortedz)
xline(78,'r--')
set(gca,'YTickLabel',{' '})
xlim([75 85])
xlabel('Validation Percentage')