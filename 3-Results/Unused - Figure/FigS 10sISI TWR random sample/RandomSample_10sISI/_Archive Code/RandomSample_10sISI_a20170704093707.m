%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GET MWT INFORMATION
N_plate = 9;
msrlist = {'RevFreq','RevSpeed','RevDur'};
complist = {'Initial', 'Hablevel', 'HabAmount'};

% get mwt paths
load(fullfile(fileparts(pSave),'MWTDB.mat'))
MWTDB_0mM = MWTDB(ismember(MWTDB.groupname,'N2'),:);
MWTDB_400mM = MWTDB(ismember(MWTDB.groupname,'N2_400mM'),:);

%% randomly select groups
A = nan(100,1);
T = table;
for msri = 1:numel(msrlist)
   T.(msrlist{msri}) = A; 
end
Initial = T;
Curve = T;
HabLevel = T;
HabPct = T;

for randsamplei = 1:3
    processIntervalReporter(3,1,'random sampling',randsamplei)
for repi = 1:100
    
    pSaveR = create_savefolder(fullfile(pSave,'sampling'),sprintf('r%ds%d',randsamplei,repi));
    
    processIntervalReporter(100,1,'-- rep',repi)
    
    % get random selected plates
    pMWT_0mM = MWTDB_0mM.mwtpath(randi(size(MWTDB_0mM,1),N_plate,1));
    pMWT_400mM = MWTDB_400mM.mwtpath(randi(size(MWTDB_400mM,1),N_plate,1));
    pMWT = [pMWT_0mM; pMWT_400mM];

    % ANALYZE DATA
    MWTSet = Dance_ShaneSpark_v1707(pMWT,pSaveR,'pctstats',false,'displayopt',false);

    %% get info
    pS = MWTSet.pSave;
    
    for msri = 1:numel(msrlist)
        msr = msrlist{msri};
        S = MWTSet.Stats.(msr);
        
        % curve
        p_curve = table2array(S.Curve.rmanova({'dose:tap'},{'pValue'}));
        Curve.(msr)(repi) = p_curve;
        
        % initial
        p_initial = S.Initial.anova.pvalue;
        p_initial_dir = S.Initial.anova.means(2) - S.Initial.anova.means(1);
        if p_initial_dir > 0; p_initial_dir = 1; else; p_initial_dir = -1; end
        Initial.(msr)(repi) = p_initial * p_initial_dir;
        
        % hab level
        p_hablevel = S.Hablevel.anova.pvalue;
        p_hablevel_dir = S.Hablevel.anova.means(2) - S.Hablevel.anova.means(1);
        if p_hablevel_dir > 0; p_hablevel_dir = 1; else; p_hablevel_dir = -1; end
        HabLevel.(msr)(repi) = p_hablevel * p_hablevel_dir;

        % hab pct
        p_hablpct = S.HabAmount.anova.pvalue;
        p_hablpct_dir = S.HabAmount.anova.means(2) - S.HabAmount.anova.means(1);
        if p_hablpct_dir > 0; p_hablpct_dir = 1; else; p_hablpct_dir = -1; end
        HabPct.(msr)(repi) = p_hablpct * p_hablpct_dir;

        
    end
end


RandomSample.Curve{randsamplei} = Curve;
RandomSample.Initial{randsamplei} = Initial;
RandomSample.HabLevel{randsamplei} = HabLevel;
RandomSample.HabPct{randsamplei} = HabPct;
end

save(fullfile(pSave,'data.mat'),'RandomSample')














