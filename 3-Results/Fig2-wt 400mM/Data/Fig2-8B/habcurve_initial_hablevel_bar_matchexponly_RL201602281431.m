%  habituation curve 10sISI of valid experiments
%% INITIALIZING
clc; clear; close all;
%% PATHS
pSave = mfilename('fullpath'); if isdir(pSave) == 0; mkdir(pSave); end
pTrvLegend = '/Users/connylin/Dropbox/RL/Code/Modules/Chor_output_handling/trv/legend_trv.mat';
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
validMWTinfo = '/Users/connylin/Dropbox/rl/Publication/PhD Dissertation/Chapters/3-STH N2/3-Results/10sISI 400mM/m1-sample/select_valid_sample_10sISI_RL201602201421/MWTInfo.mat';
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
    pmwt = fileparts(pf);
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
% measurename = {'RevDur','RevFreq','RevSpeed'};
% gname = Data.groupname;
% gu = unique(gname);
% for vi = 1:numel(measurename)
%     msr = measurename{vi};
%     for gi = 1:numel(gu)
%         fprintf('calculating %s %s\n',msr,gn);
%         gn = gu{gi};
%         i = ismember(Data.groupname,gn);
%         % calculate
%         X = Data.(msr)(i);
%         T = Data.tap(i);
%         yname = regexprep(msr,'_',' ');
%         B = table;
%         [B.tap,B.N,B.mean,B.SD,B.SE] = grpstats(X,T,{'gname','numel','mean','std','sem'});
%         B.tap = cellfun(@str2num,B.tap);
%         T = table;
%         T.groupname = repmat(gn,size(B,1),1);
%         T = [T B];
%         savename = sprintf('%s %s N by plate.csv',msr,gn);
%         cd(pSave); writetable(T,savename);
%     end
% end



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
% save
cd(pSave);
ExpMean = TS;
save('Hab exp mean.mat','ExpMean');



%% initial mean (by experiments)
measurename = {'RevDur','RevFreq','RevSpeed'};
Initial = ExpMean(ExpMean.tap == 1,:);
% include only matched experiments
a = unique(Initial.expname(ismember(Initial.groupname,'N2_400mM')));
b = unique(Initial.expname(ismember(Initial.groupname,'N2')));
c = intersect(a,b);
Initial(~ismember(Initial.expname,c),:) = [];

group = Initial.groupname;
cd(pSave);
fid = fopen('ANOVA initial.txt','w');
fprintf(fid,'ANOVA1 results for initial:\n');
for vi = 1:numel(measurename)
    msr = measurename{vi};
    fprintf(fid,'%s\n',msr);
    x = Initial.(msr);
    [p,t,s] = anova1(x,group,'off');
    result = anova_textresult(t);
    fprintf(fid,'%s\n',result);
end
fclose(fid);
% make graphing table outputs
for vi = 1:numel(measurename)
    msr = measurename{vi};
    T = descriptiveStatsTable(Initial.(msr),Initial.groupname);
    filename = sprintf('%s/Initial %s.csv',pSave,msr);
    writetable(T,filename)
end

%% initial percent control
outputname = 'Initial';
T = Initial;
measurename = {'RevDur','RevFreq','RevSpeed'};
gnu = unique(T.groupname);
enu = unique(T.expname);
fid = fopen(sprintf('%s/ttest pct %s.txt',pSave,outputname),'w');
fprintf(fid,'ttest results for %s (sig diff from 1):\n',outputname);
S = table;
S.measures = measurename'; 
S.mean = nan(size(S,1),1);
S.SE = nan(size(S,1),1);
S.N = nan(size(S,1),1);
for vi = 1:numel(measurename)
    M = table;
    M.expname = enu;
    msr = measurename{vi};
    for gi = 1:numel(gnu)
        gn= gnu{gi};
        A = T(ismember(T.groupname,gnu{gi}),:);
        [i,j] = ismember(A.expname,M.expname);
        M.(gn) = A.(msr)(j(i));
    end
    pct = M.N2_400mM./M.N2;
    [h,p,ci,s] = ttest(pct,1);
    fprintf(fid,'%s: t(%d)=%.3f, p= %.4f\n',msr,s.df,s.tstat,p);
    % stats output
    S.mean(vi) = mean(pct);
    n = numel(pct);
    S.SE(vi) = std(pct)./sqrt(n-1);
    S.N(vi) = n;
    writetable(S,sprintf('%s/%s relative2control.csv',pSave,outputname));
end
fclose(fid);

%% hab mean (by experiments)
outputname = 'HabLevel';
measurename = {'RevDur','RevFreq','RevSpeed'};
% calcualte hab level from last 3 taps
Last3Taps = ExpMean(ismember(ExpMean.tap,[28:30]),:);
group = cellfun(@strcat,Last3Taps.groupname,cellfunexpr(Last3Taps.groupname,'EXP'),Last3Taps.expname,'UniformOutput',0);
T = table;
T.gname = unique(group);
a = regexpcellout(T.gname,'EXP','split');
T.groupname = a(:,1);
T.expname = a(:,2);
for vi = 1:numel(measurename)
    msr = measurename{vi};
    [m,gn] = grpstats(Last3Taps.(msr),group,{'mean','gname'});
    [i,j] = ismember(gn,T.gname);
    T.(msr) = nan(size(T,1),1);
    T.(msr)(j(i)) = m;
end
T.gname = [];
eval(sprintf('%s=T;',outputname));

T = HabLevel;
% include only matched experiments
a = unique(T.expname(ismember(T.groupname,'N2_400mM')));
b = unique(T.expname(ismember(T.groupname,'N2')));
c = intersect(a,b);
T(~ismember(T.expname,c),:) = [];
HabLevel = T;
group = T.groupname;
cd(pSave);
fid = fopen(sprintf('ANOVA %s.txt',outputname),'w');
fprintf(fid,'ANOVA1 results for %s:\n',outputname);
for vi = 1:numel(measurename)
    msr = measurename{vi};
    fprintf(fid,'%s\n',msr);
    [p,t,s] = anova1(T.(msr),T.groupname,'off');
    result = anova_textresult(t);
    fprintf(fid,'%s\n',result);
end
fclose(fid);
% make graphing table outputs
for vi = 1:numel(measurename)
    msr = measurename{vi};
    S = descriptiveStatsTable(T.(msr),T.groupname);
    filename = sprintf('%s/%s %s.csv',pSave,msr,outputname);
    writetable(S,filename)
end

%% hab percent control
outputname = 'HabLevel';
T = HabLevel;
measurename = {'RevDur','RevFreq','RevSpeed'};
gnu = unique(T.groupname);
enu = unique(T.expname);
fid = fopen(sprintf('%s/ttest pct %s.txt',pSave,outputname),'w');
fprintf(fid,'ttest results for %s (sig diff from 1):\n',outputname);
S = table;
S.measures = measurename'; 
S.mean = nan(size(S,1),1);
S.SE = nan(size(S,1),1);
S.N = nan(size(S,1),1);
for vi = 1:numel(measurename)
    M = table;
    M.expname = enu;
    msr = measurename{vi};
    for gi = 1:numel(gnu)
        gn= gnu{gi};
        A = T(ismember(T.groupname,gnu{gi}),:);
        [i,j] = ismember(A.expname,M.expname);
        M.(gn) = A.(msr)(j(i));
    end
    pct = M.N2_400mM./M.N2;
    [h,p,ci,s] = ttest(pct,1);
    fprintf(fid,'%s: t(%d)=%.3f, p= %.4f\n',msr,s.df,s.tstat,p);
    % stats output
    S.mean(vi) = mean(pct);
    n = numel(pct);
    S.SE(vi) = std(pct)./sqrt(n-1);
    S.N(vi) = n;
        writetable(S,sprintf('%s/%s relative2control.csv',pSave,outputname));

end
fclose(fid);

















