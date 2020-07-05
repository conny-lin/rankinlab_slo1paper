%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% GLOBAL INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths ++++++++++++++
pD = fileparts(pM);
pSave = pM;
% --------------------
% load MWTDB +++++++
load(fullfile(pD,'MWTDB.mat'),'MWTDB');
pMWT = MWTDB.mwtpath;
% -----------------

% settings ++++++
% ---------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data ================================================================
% IMPORT HABITUATION CURVES DATA +++++++++++
[Data,msrlist] = importHabCurveData(pMWT);
% ------------------------------------------

% create tap data obj ++++++++++++++++++++++
TD = TapData;
v = {'mwtpath','groupname','tap'};
for vi = 1:numel(v)
    TD.(v{vi}) = Data.(v{vi});
end
TD.Response = table2array(Data(:,msrlist));
TD.Measures = msrlist;
% ------------------------------------------
% =========================================================================

%% get habituated response
tap = TD.tap;
R = array2table(TD.Response,'VariableNames',TD.Measures);
T = table;
T.mwtid = TD.mwtid;
T.tap = tap;
T = [T R];
T(~ismember(T.tap,[28:30]),:) = [];

%% summary
msrlist = TD.Measures;
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    [m,g,s] = grpstats(T.(msr),T.mwtid,{'mean','gname','sem'});
    
    T1 = table;
    T1.mwtid = cellfun(@str2num,g);
    T1.(msr) = m;
        
    if msri ==1
        T2 = T1;
    else
        T2 = outerjoin(T2,T1,'Key',{'mwtid'},'MergeKey',1);
    end
end

%% summarize table for anova
M = TD.MWTDB;
M = M(:,{'mwtid','groupname'});
gn = M.groupname;
gn = regexprep(gn,'N2_','');
gn = regexprep(gn,'mM','');
gn(ismember(gn,'N2')) = {'0'};
M.groupname = gn;
T2 = outerjoin(M,T2,'Key',{'mwtid'},'MergeKey',1);

%% anova
for msri = 1:numel(msrlist)
    msr = msrlist{msri};
    x = T2.(msr);
    group = T2.groupname;
    [text,T,p,s,t,ST] = anova1_autoresults(x,group);
    save(fullfile(pM,[msr,'.mat']),'text','T','p','s','t','ST');
    writetable(T,fullfile(pM,[msr,'.csv']));
    txt = anova1_multcomp_auto(x,group);
    fid = fopen(fullfile(pM,[msr,' ANOVA.txt']),'w');
    fprintf(fid,'%s',txt);
    fclose(fid);
end





% report done ++++++++++++
fprintf('\nDONE\n'); return
% ------------------------
% =========================================================================

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















