%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'FB','genSave',true);
% addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% load data +++++++++++++
pDataFolder = fullfile(fileparts(pM),'Dance_DrunkMoves');
datafilename = 'Dance_DrunkMoves.mat';
pData = fullfile(pDataFolder,datafilename);
load(pData,'MWTSet');
% ----------------------------

%%
DataMaster = MWTSet.Data_Plate;
msrlist = {'speed','curve'};
MWTDB = parseMWTinfo(MWTSet.PATHS.pMWT);

for msri = 1:numel(msrlist)
    % get msr name
    msr = msrlist{msri};
    % get data
    Data = DataMaster.(msr);
    m = MWTSet.Info.VarIndex.mwtname(Data.mwtname);
    g = MWTSet.Info.VarIndex.groupname(Data.groupname);
    e = MWTSet.Info.VarIndex.expname(Data.expname);
    a = cellfun(@fullfile,e,g,m,'UniformOutput',0);
    MWTDB =parseMWTinfo(a);
    MWTDBU = parseMWTinfo(unique(a));
    [i,j] = ismember(MWTDB.mwtpath, MWTDBU.mwtpath);
    Data.mwtname = j;
    
    %% modify time 
    Data.timeind(Data.timeind ==1) = 0;
    tt = [0:10:60];
    Data(~ismember(Data.timeind,tt),:) = [];
    Data.timeind = (Data.timeind./10)+1;
    
    
    %% RM ANOVA ==========================================================
    % convert data to rmANOVA format
    % get dose
    g = MWTDBU.groupname;
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
    rmTerm = sprintf('%s1-%s7~%s',rmName,rmName,gFactorName);
    rm = fitrm(T,rmTerm,'WithinDesign',withinDesignT);
    t = anovan_textresult(ranova(rm), 0, 'pvlimit',pvlimit);
    textout = sprintf('RMANOVA(%s:%s):\n%s\n',rmName,gFactorName,t);
    % ----------------------------

    % rmanova pairwise by group ++++++++
    compName = 'gname';
    rm = fitrm(T,'t1-t7~gname','WithinDesign',withinDesignT);
    t = multcompare(rm,compName);
    %  keep only unique comparisons
    a = strjoinrows([t.([compName,'_1']) t.([compName,'_2'])],' x ');
    save(fullfile(pM,[msr ,' by group.mat']),'t');
    t(~ismember(a,gpairs),:) =[];
    t(t.pValue > alpha,:) = [];
    % text output
    textout = sprintf('%s\nPosthoc(Tukey) by %s:\n',textout,compName);
    if isempty(t)
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
    if isempty(t)
        textout = sprintf('%s\nAll comparison = n.s.\n',textout);
    else
        t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
        textout = sprintf('%s\n%s\n',textout,t);
    end
    % --------------------------------
    
    % write text output ++++++++++
    fid = fopen(fullfile(pM,[msr, ' rmANOVA.txt']),'w');
    fprintf(fid,'%s',textout);
    fclose(fid);
    % ---------------------------
    %% RM ANOVA END ========================================================

    
    
    %% comparison by time within group +++++++++++
    groupnames =unique(T.gname);
    Stats = struct;
    suffix = 'btw time within group';
    fid = fopen(fullfile(pM,sprintf('%s %s.txt',msr,suffix)),'w');

    for gi = 1:numel(groupnames)
        gn = groupnames{gi};
        T1 = T(ismember(T.gname,gn),:); % get data
        % transform data
        vnames = T1.Properties.VariableNames; % get variable names
        vind = regexpcellout(vnames',rmName); % get variable names containig time
        T1 = T1(:,vind); % retain only time data
        vnames(~vind) = []; % retain only time variables
        gid = cellfun(@str2num,regexprep(vnames',rmName,''))'; % transform into numeric values
        gid(2:end) = gid(2:end)*10; % update to actual time
        rown = size(T1,1); % get row number
        gidmatrix = repmat(gid,rown,1); % make gid matrix
        gid = reshape(gidmatrix,numel(T1),1); % reshape id
        data = reshape(table2array(T1),numel(T1),1); % reshape data
        [txt,anovastats,multstats] = anova1_std_v2(data,gid);
        
        
        fprintf(fid,'*** %s ***\n',gn); % add header
        fprintf(fid,'%s\n',txt); % add to printout
    
        % collect
        Stats.(msr).(gn).anova = anovastats;
        Stats.(msr).(gn).multstats = multstats;
        return
    end
    fclose(fid); % close txt file

    % save results
    savename = sprintf('%s %s.mat',msr,suffix);
    save(fullfile(pM,savename),'Stats','MWTDB');
    % --------------------------------

    

    
end












