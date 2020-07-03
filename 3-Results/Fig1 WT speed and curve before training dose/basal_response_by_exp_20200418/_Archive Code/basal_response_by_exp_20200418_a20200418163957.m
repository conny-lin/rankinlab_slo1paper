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

%% fold 2012FP experiments into one exp
Data.expyear = regexpcellout(Data.expname,'^\d{4}','match');


%% reduce data to per plate mean
% create mwtname reference
d = Data(:,{'expyear','expname','groupname','mwtname'});
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
    

%% get N plate per exp
[t, ~, ~, l] = crosstab(Data_per_plate.expname, Data_per_plate.groupname);
T = array2table(t, ...
        'RowNames',l(1:size(t,1),1), ...
        'VariableNames', l(1:size(t,2),2));
writetable(T, fullfile(pM, 'n_plate_per_exp.csv'), 'WriteRowNames',true)

%% get N plate per exp year
[t, ~, ~, l] = crosstab(Data_per_plate.expyear, Data_per_plate.groupname);
T = array2table(t, ...
        'RowNames',l(1:size(t,1),1), ...
        'VariableNames', l(1:size(t,2),2));
writetable(T, fullfile(pM, 'n_plate_per_expyear.csv'), 'WriteRowNames',true)

%% mean per exp year per group
D = Data_per_plate;
len = size(D,1);
s = cell(len,1);
for i = 1:len
    s{i} = sprintf('%s$%s', D.expname{i}, D.groupname{i});
end
exp_gname = s;
D.exp_gname = s;

d = D(:,{'expyear','groupname','exp_gname'});
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
writetable(Data_per_exp, fullfile(pM,'desc_stats_by_expyear.csv'));


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


%% transform data to % of mean 0mM of the same experiment 
D = Data_per_exp;
i = ismember(D.groupname, 'N2');
D_ctrl = D(i,:);
D_etoh = D(~i,:);

gnameU = unique(D_etoh.groupname);

for msri = 1:numel(msrlist)
    msrname = (msrlist{msri});
    
    d = D_etoh.(msrname);
    [i,j] = ismember(D_etoh.expyear, D_ctrl.expyear);
    c = D_ctrl.(msrname)(j);
    r = d./c;
    
    D_etoh.([msrname,'_pct']) = r;
        
end

%% do stats - anova
varlist = {'speed_pct','curve_pct'};
Stats = struct;
for msri = 1:numel(varlist)
    vname = varlist{msri};
    data = D_etoh.(vname);
    group = D_etoh.groupname;
    [txt,anovastats,multstats,T,ST,ST2] = anova1_std_v2(data,group,'plate');
    
    writetable(T, fullfile(pM, sprintf('%s descriptive.csv', vname)));
    
    Stats.(vname) = ST2;
    
    fid = fopen(fullfile(pM,sprintf('%s anova.txt',vname)), 'w');
    fprintf(fid,'%s',txt);
    fclose(fid);
    
    % T test from 100% line
    fid = fopen(fullfile(pM, sprintf('%s ttest.txt',vname)),'w');
    fprintf(fid,'t test from 100%%\n');
    gu = unique(D_etoh.groupname);
    for gi = 1:numel(gu)
        gn = gu{gi};
        x = D_etoh.(vname)(ismember(D_etoh.groupname, gn));
        [h,p,c,s] = ttest(x,100);
        fprintf(fid,'%s: ttest(%d) = %.4f, p=%.4f\n',gn, s.df, s.tstat, p);
        Stats.(vname).ttest = s;
    end
    fclose(fid);
end


    

save(fullfile(pM,'data.mat'),'Data','Stats');










