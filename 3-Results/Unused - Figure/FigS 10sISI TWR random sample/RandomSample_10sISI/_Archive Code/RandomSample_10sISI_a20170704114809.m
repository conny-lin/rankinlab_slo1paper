% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% defaults
N_plate = 9;
randsampleN = 3;
repN = 100;
psig = 0.05;
msrlist = {'RevFreq','RevSpeed','RevDur'};
complist = {'Initial', 'Hablevel', 'HabAmount'};

% expected result direction
expectedDir = struct;
expectedDir.Initial.RevFreq = -1;
expectedDir.Initial.RevDur = -1;
expectedDir.Initial.RevSpeed = 1;
expectedDir.HabLevel.RevFreq = -1;
expectedDir.HabLevel.RevDur = 1;
expectedDir.HabLevel.RevSpeed = 1;
expectedDir.HabPct.RevFreq = 1;
expectedDir.HabPct.RevDur = -1;
expectedDir.HabPct.RevSpeed = -1;

% get mwt paths
load(fullfile(fileparts(pSave),'MWTDB.mat'))
MWTDB_0mM = MWTDB(ismember(MWTDB.groupname,'N2'),:);
MWTDB_400mM = MWTDB(ismember(MWTDB.groupname,'N2_400mM'),:);

% create template output
A = nan(100,1);
T = table;
for msri = 1:numel(msrlist)
   T.(msrlist{msri}) = A; 
end
Initial = T;
Curve = T;
HabLevel = T;
HabPct = T;

%% randomly select groups
for randsamplei = 1:randsampleN
    
    for repi = 1:repN

        pSaveR = create_savefolder(fullfile(pSave,'sampling'),sprintf('r%ds%d',randsamplei,repi));
        processIntervalReporter(100,1,sprintf('-- sample%d rep',randsamplei),repi);

        % get random selected plates
        pMWT_0mM = MWTDB_0mM.mwtpath(randi(size(MWTDB_0mM,1),N_plate,1));
        pMWT_400mM = MWTDB_400mM.mwtpath(randi(size(MWTDB_400mM,1),N_plate,1));
        pMWT = [pMWT_0mM; pMWT_400mM];

        % ANALYZE DATA
        if ~isdir(fullfile(pSaveR,'Dance_ShaneSpark_v1707'))
            MWTSet = Dance_ShaneSpark_v1707(pMWT,pSaveR,'pctstats',false,'displayopt',false);
        else
            load(fullfile(pSaveR,'Dance_ShaneSpark_v1707','Dance_ShaneSpark_v1707.mat'));
        end

        % get info
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
    
    % store in struct
    RandomSample.Curve{randsamplei} = Curve;
    RandomSample.Initial{randsamplei} = Initial;
    RandomSample.HabLevel{randsamplei} = HabLevel;
    RandomSample.HabPct{randsamplei} = HabPct;
end
save(fullfile(pSave,'data.mat'),'RandomSample');


%% create summary table
component_names = fieldnames(RandomSample);
T = table;
for randsamplei = 1:randsampleN
    L = {};
    Array = {};
    for cmpnamei = 2:numel(component_names)
        cmpname = component_names{cmpnamei};
        D = table2array(RandomSample.(cmpname){randsamplei});

        if ~strcmp(cmpname,'Curve')            
            % get expected value
            eValue = nan(size(msrlist));
            for msri = 1:numel(msrlist)
                eValue(msri) = expectedDir.(cmpname).(msrlist{msri});
            end
            D = repmat(eValue,size(D,1),1).*D;
        end
        % get percent pass pvalue
        P = D<psig; % get pvalue valid
        ps = sum(P)./repN; % get percentage
        % construct output
        m = msrlist';
        c = repmat({cmpname},numel(m),1);
        p = ps';
        
        Array{cmpnamei,1} = p;
        L = [L;[m c]];
    end
    Array = cell2mat(Array);
    vname = sprintf('T%d',randsamplei);
    if randsamplei ==1
        T = cell2table(L,'variableNames',{'msr','Components'});
        T.(vname) = Array;
    else
        T1 = cell2table(L,'variableNames',{'msr','Components'});
        T1.(vname) = Array;
        T = innerjoin(T, T1);
    end
end

writetable(T,fullfile(pSave,'random sample result.csv'));


















