%  habituation curve 10sISI of valid experiments

%% INITIALIZING
clc; clear; close all;
%% PATHS
validMWTinfo = '/Users/connylin/Dropbox/RL/Publication InPrep/PhD Dissertation/Chapters/3-STH N2/2-Method/3-valid samples/select_valid_sample_10sISI_RL201602201421/MWTInfo.mat';
pSave = mfilename('fullpath'); if isdir(pSave) == 0; mkdir(pSave); end
%% ADD FUNCTION PATH
pFun = {'/Users/connylin/Dropbox/Code/Matlab/Library/General';
        '/Users/connylin/Dropbox/RL/Code/Modules/MWTDatabase';
        '/Users/connylin/Dropbox/Code/Matlab/Library/Graphs';
        '/Users/connylin/Dropbox/Code/Matlab/Library/Stats';
        '/Users/connylin/Dropbox/RL/Code/Modules/Chor'};
for x = 1:numel(pFun); 
    p = pFun{x};
    if isdir(p) == 0; error('function path %d not found'); end
    addpath(p); 
end
%% LOAD MWT DATABASE
load(validMWTinfo);
pMWT = MWTDB.mwtpath;
%% CHOR FILES
% check trv files
[pTrv,pMWThasTRV,pMWTnoTRV] = getpath2chorfile(pMWT,'*.trv');
% chor
[Legend,pMWTcS,fval,chorscript] = chormaster4('BeethovenOnly',pMWTnoTRV);
% check version of trv files
trv_oldversion = true(size(pTrv));
for mwti =1:numel(pTrv)
    pf = pTrv{mwti};
    % see version of trv
    fileID = fopen(pf,'r');
    a = textscan(fileID,'%s', 1,'Delimiter', '', 'WhiteSpace', '');
    a = char(a{1}(1));
    fc = a(1);
    fclose(fileID);
    if regexp(fc,'\d{1,}') == 1; 
        trv_oldversion(mwti) = false;
    end
end
fprintf('%d files are old trv version\n',sum(trv_oldversion))




%% HABITUATION CURVES
% check trv files
[pTrv,pMWThasTRV,pMWTnoTRV] = getpath2chorfile(pMWT,'*.trv');
% import legend
pTrvLegend = '/Users/connylin/Dropbox/RL/Code/Modules/Chor_output_handling/trv/legend_trv.mat';
load(pTrvLegend);
disp(char(legend_trv'));
% get necessary legends
legend_output = {'tap','time', 'N_alreadyRev', 'N_ForwardOrPause', 'N_Rev', 'RevDis', 'RevDur'};
ind_get = find(ismember(legend_trv,legend_output));

%% import into structural array
% Import= cell(size(pTrv));
S = struct;
for mwti =1:numel(pTrv)
    pf = pTrv{mwti};
    if size(pf,1) > 1; error('fix trv paths'); end
    % see version of trv
    d = dlmread(pf);
    ncol = size(d,2);
    if ncol ~= numel(legend_trv); error('trv col number wrong'); end
    % add tap
    tapnumber = (1:size(d,1))';
    D = array2table([tapnumber,d(:,ind_get)],'VariableNames',legend_output);
    % frequency is # reversed divided by total number not already reversing
    D.RevFreq = D.N_Rev./(D.N_Rev + D.N_ForwardOrPause);
    % calculate reversal speed: reversal speed is reversal distance divided by reversal duration
    D.RevSpeed = D.RevDis./D.RevDur;
    D.RevDis = [];
    % prepare summary output identifiers
    db = parseMWTinfo(fileparts(pf));
    dbname = {'expname','groupname','mwtname'};
    S(mwti).mwtpath = {pmwt};
    S(mwti).expname = db.expname;
    S(mwti).groupname = db.groupname;
    S(mwti).mwtname =  db.mwtname;
    % enter variable in structural array
    vnames = D.Properties.VariableNames;
    for vi = 1:numel(vnames)
        vn = vnames{vi};
        S(mwti).(vn) = D.(vn);
    end
end

%% expand import into table
n = arrayfun(@(x) numel(x.tap),S);
nrow = sum(n);
rowi = [0 cumsum(n)];
vnames = fieldnames(S);
identifiers = {'expname', 'groupname', 'mwtname','mwtpath'};
varn = vnames(~ismember(vnames,identifiers));
measurename = varn(~ismember(varn,'wormid'));
% prepare table
Data = table;
for vi = 1:numel(identifiers)
    v = identifiers{vi};
    a = cell(nrow,1); 
    for i = 1:numel(S); a(rowi(i)+1:rowi(i+1)) = S(i).(v); end
    Data.(v) = a;
end
% variables
for vi = 1:numel(varn)
    v = varn{vi};
    A  = cell(size(S)); for i = 1:numel(S); A{i} = S(i).(v); end
    Data.(v) = cell2mat(A');
end


%% create stats output for hab curve by plates
measurename = {'RevDur','RevFreq','RevSpeed'};
gname = Data.groupname;
gu = unique(gname);
for vi = 1:numel(measurename)
    msr = measurename{vi};
    for gi = 1:numel(gu)
        fprintf('calculating %s %s\n',msr,gn);
        gn = gu{gi};
        i = ismember(Data.groupname,gn);
        % calculate
        X = Data.(msr)(i);
        T = Data.tap(i);
        yname = regexprep(msr,'_',' ');
        B = table;
        [B.tap,B.N,B.mean,B.SD,B.SE] = grpstats(X,T,{'gname','numel','mean','std','sem'});
        B.tap = cellfun(@str2num,B.tap);
        T = table;
        T.groupname = repmat(gn,size(B,1),1);
        T = [T B];
        savename = sprintf('%s %s N by plate.csv',msr,gn);
        cd(pSave); writetable(T,savename);
    end
end


%% create stats output for hab curve by experiments
% calculate experiment mean
measurename = {'RevDur','RevFreq','RevSpeed'};
gname = Data.groupname;
gu = unique(gname);
% convert group name into label name
a = Data.groupname; a(ismember(a,'N2')) = {'0mM'}; a = regexprep(a,'(N2_)|(mM)','');
gnamelabel = a;

for vi = 1:numel(measurename)
    msr = measurename{vi};
    T1 = table;
    for gi = 1:numel(gu)
        gn = gu{gi};

        fprintf('calculating %s %s\n',msr,gn);
        i = ismember(Data.groupname,gn);
        eu = unique(Data.expname(i));
        for ei = 1:numel(eu)
            en = eu{ei};
            i = ismember(Data.groupname,gn) & ismember(Data.expname,en);
            % calculate
            X = Data.(msr)(i);
            tap = Data.tap(i);
            B = table;
            [B.tap,B.N,B.(msr)] = grpstats(X,tap,{'gname','numel','mean'});
            B.tap = cellfun(@str2num,B.tap);
            T = table;
            T.groupname = repmat({gn},size(B,1),1);
            T.expname = repmat({en},size(B,1),1);
            T = [T B];
            T1 = [T1;T];
        end
    end
    if vi ==1;
        TS = T1;
    else
        TS = innerjoin(TS,T1);
    end
end

%% output by experiments
measurename = {'RevDur','RevFreq','RevSpeed'};
gname = TS.groupname;
gu = unique(gname);
for vi = 1:numel(measurename)
    msr = measurename{vi};
    for gi = 1:numel(gu)
        fprintf('calculating %s %s\n',msr,gn);
        gn = gu{gi};
        i = ismember(TS.groupname,gn);
        % calculate
        X = TS.(msr)(i);
        T = TS.tap(i);
        yname = regexprep(msr,'_',' ');
        B = table;
        [B.tap,B.N,B.mean,B.SD,B.SE] = grpstats(X,T,{'gname','numel','mean','std','sem'});
        B.tap = cellfun(@str2num,B.tap);
        T = table;
        T.groupname = repmat(gn,size(B,1),1);
        T = [T B];
        savename = sprintf('%s %s N by experiment.csv',msr,gn);
        cd(pSave); writetable(T,savename);
    end
end


%% repeated measures ANOVA by experiment
measurename = {'RevDur','RevFreq','RevSpeed'};
gname = TS.groupname;
gu = unique(gname);
for vi = 1:numel(measurename)
    msr = measurename{vi};
    X = cell(numel(gu));
    for gi = 1:numel(gu)
        fprintf('calculating %s %s\n',msr,gn);
        gn = gu{gi};
        i = ismember(TS.groupname,gn);
        % calculate
        x = TS.(msr)(i);
        tap = TS.tap(i);
        en = TS.expname(i);
        ntap = numel(unique(tap));
        enu = unique(en);
        nexp = numel(enu);
        a = nan(nexp,ntap);
        for ei = 1:numel(enu)
            e = enu(ei);
            i = ismember(en,e);
            a(ei,tap(i)) = x(i)';
        end
        i = ~any(isnan(a),2);
        n = sum(i);
        fprintf('N exp = %d\n',n);
        X{gi} = a(i,:);
    end
    [p,t]  = anova_rm(X,'off');
    a = t(2:end,2:end);
    colname =  regexprep(t(1,2:end),'>F','');
    rowname = regexprep(t(2:end,1),'[(]|[)]','');
    a = array2table(a,'VariableNames',colname);
    A = table;
    A.rowname = rowname;
    A = [A a];
    cd(pSave);
    writetable(A,sprintf('%s RMANOVA.csv',msr))
end




















