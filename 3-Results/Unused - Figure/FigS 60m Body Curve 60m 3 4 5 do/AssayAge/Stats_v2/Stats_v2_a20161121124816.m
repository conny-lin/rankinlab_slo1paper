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
DataMaster = MWTSet.Data_Plate;
msrlist = {'speed','curve'};




%%

MWTDB = parseMWTinfo(MWTSet.PATHS.pMWT);

for msri = 1:numel(msrlist)
    % get msr name
    msr = msrlist{msri};
    % get data
    Data = DataMaster.(msr);
    Data.groupname = MWTSet.Info.VarIndex.groupname(Data.groupname);
    Data.expname = MWTSet.Info.VarIndex.expname(Data.expname);

    %% RM ANOVA ==========================================================
    % convert data to rmANOVA format
    % get dose
    g = Data.groupname;
    a = regexpcellout(g,'\d{1,}(?=mM)','match');
    a(cellfun(@isempty,a)) = {'0'};
    dose = cellfun(@str2num,a);
    
    % get age
    a = regexpcellout(g,'\d{1,}(?=d)','match');
    a(cellfun(@isempty,a)) = {'4'};
    age = cellfun(@str2num,a);
    
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
    
    
    %% combine all
    T = table;
    T.gname = g;
    T.age = age;
    T.dose = dose;
    T = [T t];
    
    % exclude no data or zero data
    T(any(isnan(A),2),:) = [];
    T(any(A==0,2),:) = [];
    
    % within design table
    withinDesignT = table;
    withinDesignT.t = unique(Data.timeind);
    
    % get unique pairwise
    gpairs = pairwisecomp_getpairs(unique(T.gname));
    gpairs = strjoinrows(gpairs,' x ');

    % rmanova settigns +++
    pvlimit = 0.001;
    alpha = 0.05;
    rmName = 't';    
    gFactorName = 'dose*age';
    % ------------------
    
    % rmanova multi-factor +++++++
    rmTerm = sprintf('%s1-%s60~%s',rmName,rmName,gFactorName);
    rm = fitrm(T,rmTerm);
    t = anovan_textresult(ranova(rm), 0, 'pvlimit',pvlimit);
    textout = sprintf('RMANOVA(%s:%s):\n%s\n',rmName,gFactorName,t);
    % ----------------------------

    % rmanova pairwise by group ++++++++
    compName = 'gname';
    rm = fitrm(T,'t1-t60~gname','WithinDesign',withinDesignT);
    t = multcompare(rm,compName);
    %  keep only unique comparisons
    a = strjoinrows([t.([compName,'_1']) t.([compName,'_2'])],' x ');
    save(fullfile(pM,[msr ,' by group.mat']),'t');
    t(~ismember(a,gpairs),:) =[];
    t(t.pValue > alpha,:) = [];
    % text output
    textout = sprintf('%s\nPosthoc(Tukey) by %s:\n',textout,compName);
    if isempty(t); 
        textout = sprintf('%s\nAll comparison = n.s.\n',textout);
    else
        t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
        textout = sprintf('%s\n%s\n',textout,t);
    end
    % ----------------------------------

    % comparison by time +++++++++++
    compName = 'gname';
    t = multcompare(rm,compName,'By',rmName);
    
    %  keep only unique comparisons
    a = strjoinrows([t.([compName,'_1']) t.([compName,'_2'])],' x ');
    save(fullfile(pM,[msr ,' by time.mat']),'t');
    t(~ismember(a,gpairs),:) =[];
    t(t.pValue > alpha,:) = [];
    % record
    textout = sprintf('%s\n\nPosthoc(Tukey)%s by %s:\n',textout,rmName,compName); 
    if isempty(t); 
        textout = sprintf('%s\nAll comparison = n.s.\n',textout);
    else
        t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
        textout = sprintf('%s\n%s\n',textout,t);
    end
    % --------------------------------
    %% RM ANOVA END ========================================================

    
    %% write text output ++++++++++
    fid = fopen(fullfile(pM,[msr, ' rmANOVA.txt']),'w');
    fprintf(fid,'%s',textout);
    fclose(fid);
    %% ---------------------------
    
end












