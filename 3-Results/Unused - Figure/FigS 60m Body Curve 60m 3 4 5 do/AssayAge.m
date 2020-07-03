%% OBJECTIVE:
% - create a list of experiments relevant for rapid tolerance chapter
%% INITIALIZING
clc; clear; close all;
%% PATHS
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);

%% load database
pDB = '/Volumes/COBOLT/MWT/MWTDB.mat';
load(pDB);
MWTDB = MWTDB.text;

%% target list
rclist = {'3600s0x0s0s'};

%% query
i = ismember(MWTDB.rc, rclist) & ismember(MWTDB.strain,'N2')...
    & ~ismember(MWTDB.groupname, {'N2_Test','N2_400mM_NoFood','N2_NoFood'});
MWTDB(~i,:) = [];
%% change to group name alternate
cd('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/MWTDatabase')
rxname = readtable('rxname_alternate.csv');
for ri =1:size(rxname,1)
    i= ismember(MWTDB.rx,rxname.rx(ri));
    MWTDB.rx(i) = rxname.rx_alternate(ri);
end
% create new names
A = [MWTDB.strain MWTDB.rx];
MWTDB.groupname = strjoinrows(A,'_');

%% keep only 1 hour recovery time
% Rx = parseRx(MWTDB.rx);
% MWTDB(Rx.rec_hr~=1,:) = [];
unique(MWTDB.groupname)

%% run Dance
pMWT = MWTDB.mwtpath;
MWTSet = Dance_DrunkMoves(pMWT,'pSave',pSave);

%% get SE
msrlist = fieldnames(MWTSet.Data_GroupByPlate);
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    T = MWTSet.Data_GroupByPlate.(msr);
    T.groupname = MWTSet.Info.VarIndex.groupname(T.groupname);
    
    gnu = unique(T.groupname);
    tn = unique(T.timeind);
    t = nan(numel(gnu),numel(tn));
    for gi = 1:numel(gnu)
        i = ismember(T.groupname,gnu(gi));
        t(gi,:) = T.SE(i)';
    end
    t = array2table(t);
    T1 = table;
    T1.groupname = gnu;
    T1 = [T1 t];
    cd(pSave);
    writetable(T1,sprintf('%s/Dance_DrunkMoves/Data Curve/%s SE.csv',pSave,msr));
end

%% run stats on 35min mark, speed and curve
pSaveA = sprintf('%s/Stats',pSave);
if isdir(pSaveA) == 0; mkdir(pSaveA); end
gnind = MWTSet.Info.VarIndex.groupname;
strainind = MWTSet.Info.VarIndex.strain;
msrlist = {'speed','curve'};
timeT = 35;
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    D = MWTSet.Data_Plate.(msr);
    D(D.timeind~=timeT,:) = [];
    gn = gnind(D.groupname);
    % age
    a = regexpcellout(gn,'\d{1,}(?=d)','match');
    a(cellfun(@isempty,a)) = {'4'};
    age = a;
    % dose
    a = regexpcellout(gn,'\d{1,}(?=mM)','match');
    a(cellfun(@isempty,a)) = {'0'};
    dose = a;
    %% anovav
    Y = D.mean;
    G = [{age} {dose}];
    anovan_std(Y,G,{'age','dose'},pSaveA,msr)    
end



%% percent diff from 0mM control
%% create exp mean
MWTSet.Data_Exp = struct;
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    D = MWTSet.Data_Plate.(msr);
    enu = unique(D.expname);
    T = table;
    for ei = 1:numel(enu)
        D1 = D(D.expname==enu(ei),:);
        [m,gn] = grpstats(D1.mean,D1.groupname,{'mean','gname'});
        T1 = table;
        T1.expname = repmat(enu(ei),numel(m),1);
        T1.groupname = cellfun(@str2num,gn);
        T1.mean = m;
        T = [T ;T1];
        
    end
    MWTSet.Data_Exp.(msr) = T;
end

%% get percent diff
gnind = MWTSet.Info.VarIndex.groupname;
strainind = MWTSet.Info.VarIndex.strain;
msrlist = {'speed','curve'};
timeT = 35;
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    D = MWTSet.Data_Plate.(msr);
    D(D.timeind~=timeT,:) = [];
    gn = gnind(D.groupname);
    % age
    a = regexpcellout(gn,'\d{1,}(?=d)','match');
    a(cellfun(@isempty,a)) = {'4'};
    age = a;
    % dose
    a = regexpcellout(gn,'\d{1,}(?=mM)','match');
    a(cellfun(@isempty,a)) = {'0'};
    dose = a;
    % get control name
    gex = regexpcellout(gn,'[_]\d{1,}mM','match');
    gnc = cell(size(gex));
    for gi = 1:numel(gex)
       gnc(gi) = regexprep(gn(gi),gex(gi),'');
    end
    [i,gnci] = ismember(gnc,gnind);
    % get control values
    Dref = MWTSet.Data_Exp.(msr);
    D.ctrl = nan(size(D,1),1);
    for di = 1:size(D,1)
        d = D(di,:);
        i = Dref.expname == d.expname & Dref.groupname == gnci(di);
        ind = find(i);
        if numel(ind)>1; error('wrong'); end
        D.ctrl(di) = Dref.mean(ind);
    end
    % calculate percent diff
    D.pct = pct_diff(D.mean, D.ctrl);
    
    % export graphing values
    T = grpstatsTable(D.pct,gnind(D.groupname));
    writetable(T,sprintf('%s/%s pct.csv',pSaveA,msr));
    
    % anovav
    Y = D.pct;
    G = [{age} {dose}];
    anovan_std(Y,G,{'age','dose'},pSaveA,msr,'suffix','pct')    
end








































