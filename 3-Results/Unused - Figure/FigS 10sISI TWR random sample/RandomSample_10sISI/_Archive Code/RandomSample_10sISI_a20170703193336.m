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

for repi = 1:100
    
    pSaveR = create_savefolder(pSave,['r',num2str(repi)]);
    
    % get random selected plates
    pMWT_0mM = MWTDB_0mM.mwtpath(randi(size(MWTDB_0mM,1),N_plate,1));
    pMWT_400mM = MWTDB_400mM.mwtpath(randi(size(MWTDB_400mM,1),N_plate,1));
    pMWT = [pMWT_0mM; pMWT_400mM];

    % ANALYZE DATA
    MWTSet = Dance_ShaneSpark_v5r1(pMWT,'pSave',pSaveR,'pctstats',false);

    %% get info
    pS = MWTSet.PATHS.pSaveA;
    
    for msri = 1:numel(msrlist)
        msr = msrlist{msri};
        S = MWTSet.Stats.(msr);
        
        % curve
        fid = fopen(fullfile(pS,'RMANOVA.txt'),'r');
        d = textscan(fid, '%s%[^\n\r]', 'Delimiter', '',  'ReturnOnError', false);
        fclose(fid);
        d = d{1};
        sec = find(regexpcellout(d,'^[-]{5}'));
        headeri = find(regexpcellout(d,sprintf('----- %s -----',msr)));
        k = sec(find((sec-headeri)==0));
        if (k+1) > numel(sec)
            j = numel(d);
        else
            j = sec(k+1)-1;
        end
        d = d(k:j);

        a = d((regexpcellout(d,'^F[(]\d{1,},\d{1,}[)]dose[*]tap')));
        % a = regexpcellout(a,', ','split');
        p = regexpcellout(a,'(?<=p (=|<))\d{1,}','match');
        p = p{1};
        if isempty(p)
            p=1;
        else
            p=str2double(p);
        end
        p_curve = p;
        Curve.(msr)(repi) = p_curve;
        
        % initial
        p_initial = S.Initial.anova.pvalue;
        Initial.(msr)(repi) = p_initial;
        
        % hab level
        p_hablevel = S.Hablevel.anova.pvalue;
        HabLevel.(msr)(repi) = p_hablevel;

        % hab pct
        p_hablpct = S.HabAmount.anova.pvalue;
        HabPct.(msr)(repi) = p_hablpct;

        
    end
end

return
















