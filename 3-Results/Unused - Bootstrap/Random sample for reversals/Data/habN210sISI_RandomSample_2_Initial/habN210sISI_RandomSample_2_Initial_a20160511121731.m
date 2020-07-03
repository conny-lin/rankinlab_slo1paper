%  habituation curve 10sISI of valid experiments

%% INITIALIZING
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
pSave = fileparts(pM);

%% LOAD MWT DATABASE
load([pM,'/DataTrv.mat']);

return


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




















