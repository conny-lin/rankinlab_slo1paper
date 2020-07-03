%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'FB','genSave',true);
% addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% load data +++++++++++++
pData = '/Users/connylin/Dropbox/RL Pub PhD Dissertation/Chapters/2-Wildtype/3-Results/0-Data/0min exposure 60min/Age N2 3 4 5 comparison/AssayAge/Dance_DrunkMoves/Dance_DrunkMoves.mat';
load(pData);
% ----------------------------

%%
DataMater = MWTSet.Data_Plate;
msrlist = {'speed','curve'};
MWTDB = parseMWTinfo(MWTSet.PATHS.pMWT);

for msri = 1:numel(msrlist)
    % get msr name
    msr = msrlist{msri};
    % get data
    Data = DataMater.(msr);
    
    %% RM ANOVA ==========================================================
    %% convert data to rmANOVA format
    % make gfactor ++++++++++++++
    gn = MWTDB.groupname;
    groupnames = gn(Data.mwtname);
    % get dose
    a = regexpcellout(groupnames,'\d{1,}(?=mM)','match');
    a(cellfun(@isempty,a)) = {'0'};
    dose = a;
    % get age
    a = regexpcellout(groupnames,'\d{1,}(?=d)','match');
    a(cellfun(@isempty,a)) = {'4'};
    age = a;
%     % make entry table
%     gfactors = table;
%     gfactors.dose = dose;
%     gfactors.age = age;
%     % ---------------------------
%     
%     % make rm factor +++++++++
%     rmfactors = table;
%     rmfactors.t = Data.timeind;
%     % ---------------------------
    
    %% table
    DataT = table;
    DataT.dose = dose;
    DataT.age = age;
    DataT.t = Data.timeind;
    
    withinDesignT = table;
    withinDesignT.t = unique(Data.t);
    

    %%
    rm = fitrm(DataT,rmTerms,'WithinDesign',tptable);

    return
    
    
%     textout = anovarm_std(T,'t',{'age','dose'},'groupnames')
    
    % rmanova multi-factor +++++++
    rmTerms = sprintf('%s%d-%s%d~%s',rmName,tpu(1),rmName,tpu(end),factorName);
    rm = fitrm(D,rmTerms,'WithinDesign',tptable);
    t = anovan_textresult(ranova(rm), 0, 'pvlimit',pvlimit);
    textout = sprintf('%s\nRMANOVA(%s:%s):\n%s\n',textout,rmName,factorName,t);
    % ----------------------------

    % rmanova pairwise ++++++++
    rmTerms = sprintf('%s%d-%s%d~%s',rmName,tpu(1),rmName,tpu(end),compName);
    rm = fitrm(D,rmTerms,'WithinDesign',tptable);
    t = multcompare(rm,compName);
    % text output
    textout = sprintf('%s\nPosthoc(Tukey) by %s:\n',textout,compName);
    if isempty(t); 
        textout = sprintf('%s\nAll comparison = n.s.\n',textout);
    else
        t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
        textout = sprintf('%s\n%s\n',textout,t);
    end
    % ----------------------------------

    % comparison by taps +++++++++++
    t = multcompare(rm,compName,'By',rmName);
    %  keep only unique comparisons
    a = strjoinrows([t.([compName,'_1']) t.([compName,'_2'])],' x ');
    t(~ismember(a,gpairs),:) =[]; 
    % record
        textout = sprintf('%s\n\nPosthoc(Tukey)%s by %s:\n',textout,rmName,compName); 
    if isempty(t); 
        textout = sprintf('%s\nAll comparison = n.s.\n',textout);
    else
        t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
        textout = sprintf('%s\n%s\n',textout,t);
    end
    % --------------------------------
    
    return
end