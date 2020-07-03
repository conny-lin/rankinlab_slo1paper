%% To analyze body curves % difference N = experiments

%% INITIALIZING -----------------------------------------------------------
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
% addpath(pM);
% -------------------------------------------------------------------------


%% load data --------------------------------------------------------------
pDataFolder = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Fig1 WT speed and curve before training dose/BasalResponse_Dose_RL201602210653';
datafilename = 'data.mat';
pData = fullfile(pDataFolder,datafilename);
load(pData);
% -------------------------------------------------------------------------


%% global variables -------------------------------------------------------
msrlist = {'speed','curve'};
% -------------------------------------------------------------------------


%% reduce data to per plate mean

% create mwtname reference
d = Data(:,{'expname','groupname','mwtname'});
mwtname_attributes = unique(d, 'rows');

T = table;
for msrlisti = 1:numel(msrlist)
    msrname = msrlist{msrlisti};
    x = Data.(msrname);
    g = Data.mwtname;
    T1 = table;

    [T1.mwtname, T1.wormN, T1.(msrname), T1.([msrname,'_SE'])] = ...
          grpstats(x,g,{'gname','numel','mean','sem'});
    if msrlisti > 1
        T = outerjoin(T,T1,'Keys',{'mwtname','wormN'}, 'MergeKeys', true);
    else
        T = T1;
        
    end
end

Data_per_plate = innerjoin(mwtname_attributes, T, 'Keys', {'mwtname'});

writetable(Data_per_plate, fullfile(pM,'desc_stats_by_plate.csv'))
    

%% mean per exp per group
D = Data_per_plate;
len = size(D,1);
s = cell(len,1);
for i = 1:len
    s{i} = sprintf('%s$%s', D.expname{i}, D.groupname{i});
end
exp_gname = s;
D.exp_gname = s;

d = D(:,{'expname','groupname','exp_gname'});
var_attribues = unique(d, 'rows');

T = table;
for msrlisti = 1:numel(msrlist)
    msrname = msrlist{msrlisti};
    x = D.(msrname);
    g = exp_gname;
    
    T1 = table;
    [T1.exp_gname, T1.plateN, T1.(msrname), T1.([msrname,'_SE'])] = ...
          grpstats(x,g,{'gname','numel','mean','sem'});
     
    if msrlisti > 1
        T = outerjoin(T,T1,'Keys',{'exp_gname','plateN'}, 'MergeKeys', true);
    else
        T = T1;
        
    end
end

Data_per_exp = innerjoin(var_attribues, T, 'Keys', {'exp_gname'});

writetable(Data_per_exp, fullfile(pM,'desc_stats_by_exp.csv'));



%% get mean per group
D = Data_per_exp;

d = D(:,{'groupname'});
var_attribues = unique(d, 'rows');

T = table;
for msrlisti = 1:numel(msrlist)
    msrname = msrlist{msrlisti};
    x = D.(msrname);
    g = D.groupname;
    
    T1 = table;
    [T1.groupname, T1.plateN, T1.(msrname), T1.([msrname,'_SE'])] = ...
          grpstats(x,g,{'gname','numel','mean','sem'});
     
    if msrlisti > 1
        T = outerjoin(T,T1,'Keys',{'groupname','plateN'}, 'MergeKeys', true);
    else
        T = T1;
        
    end
end

Data_per_group = innerjoin(var_attribues, T, 'Keys', {'groupname'});

writetable(Data_per_group, fullfile(pM,'desc_stats_by_group.csv'));


return

%% transform data to % of mean 0mM of the same experiment 


%% transform % change data to per exp

%% do stats


%% settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msrlist = {'speed','curve'};
incr = 5; % increment
tt = [5:incr:60]; % decide time frames to assay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MWTDB = parseMWTinfo(MWTSet.PATHS.pMWT);
for msri = 1:numel(msrlist)
    % get msr name
    msr = msrlist{msri};
    
    data = DataMaster.(msr);
    groups = DataMaster.groupname;
    [txt,anovastats,multstats] = anova1_std_v2(data,groups)
    
    return

    
    %% comparison by time within group +++++++++++
    groupnames =unique(T.gname);
    Stats = struct;
    suffix = 'btw time within group';
    fid = fopen(fullfile(pM,sprintf('%s %s.txt',msr,suffix)),'w');
    fid2 = fopen(fullfile(pM,sprintf('%s %s neighboring time only.txt',msr,suffix)),'w');

    for gi = 1:numel(groupnames)
        gn = groupnames{gi};
        T1 = T(ismember(T.gname,gn),:); % get data
        % transform data
        vnames = T1.Properties.VariableNames; % get variable names
        vind = regexpcellout(vnames',rmName); % get variable names containig time
        T1 = T1(:,vind); % retain only time data
        vnames(~vind) = []; % retain only time variables
        gid = cellfun(@str2num,regexprep(vnames',rmName,''))'; % transform into numeric values
        gid = gid*incr; % update to actual time
        rown = size(T1,1); % get row number
        gidmatrix = repmat(gid,rown,1); % make gid matrix
        gid = reshape(gidmatrix,numel(T1),1); % reshape id
        data = reshape(table2array(T1),numel(T1),1); % reshape data
        [txt,anovastats,multstats] = anova1_std_v2(data,gid);
        % central txt output
        fprintf(fid,'*** %s ***\n',gn); % add header
        fprintf(fid,'%s\n',txt); % add to printout
        
        % excel matrix output
        M = nan(timeN);
        r = cellfun(@str2num,multstats.gname_1);
        c = cellfun(@str2num,multstats.gname_2);
        ttr = tt; 
        [i,jr] = ismember(r,ttr);
        [i,jc] = ismember(c,ttr);
        idx = sub2ind(size(M),jr,jc);
        M(idx) = multstats.pValue;
        MTable = array2table(M);
        leg = num2cellstr(ttr');
        MTable.Properties.RowNames = leg; % create row names
        MTable.Properties.VariableNames = strjoinrows([repmat({'t'},numel(leg),1) leg],''); % create col names
        fname = fullfile(pM,sprintf('%s btw time within group pvalue Matrix %s.csv',msr, gn));
        writetable(MTable,fname,'WriteRowNames',1);
        
        % find maxed out change stats
        a = strjoinrows([num2cellstr(ttr(1:end-1))' num2cellstr(ttr(2:end))'],'*'); % create list of times
        b = strjoinrows([multstats.gname_1,multstats.gname_2],'*'); % create list of times from stats results
        A = multstats(ismember(b,a),:); % get neighboring data
        r = multcompare_text2016b(A); 
        % txt output
        fprintf(fid2,'*** %s ***\n',gn); % add header
        fprintf(fid2,'%s\n',r); % add to printout
        
        
        % collect to data output
        Stats.(msr).(gn).anova = anovastats;
        Stats.(msr).(gn).multstats = multstats;
        
    end
    fclose(fid); % close txt file

    % save results
    savename = sprintf('%s %s.mat',msr,suffix);
    save(fullfile(pM,savename),'Stats','MWTDB');
    % --------------------------------

    

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











