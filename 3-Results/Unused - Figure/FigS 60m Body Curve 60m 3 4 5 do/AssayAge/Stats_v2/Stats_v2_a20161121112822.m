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
    % get dose
    g = MWTDB.groupname;
    a = regexpcellout(g,'\d{1,}(?=mM)','match');
    a(cellfun(@isempty,a)) = {'0'};
    dose = a;
    
    % get age
    a = regexpcellout(g,'\d{1,}(?=d)','match');
    a(cellfun(@isempty,a)) = {'4'};
    age = a;
    
    % get rm variables
    mwtid = unique(Data.mwtname);
    timeid = unique(Data.timeind);
    
    A = nan(max(mwtid),max(timeid));
    i = sub2ind(size(A), Data.mwtname, Data.timeind);
    A(i) = Data.mean;
    a = cellfun(@num2str,num2cell(timeid),'UniformOutput',0);
    b = cellfunexpr(a,'t');
    c = strjoinrows([b a],'');
    t = array2table(A,'VariableNames',c);
    
    
    % combine all
    T = table;
    T.gname = g;
    T.age = age;
    T.dose = dose;
    T = [T t];
    T(any(isnan(A),2),:) = [];
    
    % within design table
    withinDesignT = table;
    withinDesignT.t = unique(Data.timeind);
    

    % run rm anova
    rm = fitrm(T,'t1-t60~dose*age','WithinDesign',withinDesignT)

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