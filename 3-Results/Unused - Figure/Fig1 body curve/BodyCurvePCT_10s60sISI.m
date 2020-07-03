%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170726
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170726



%% GET DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170726
pData10 = '/Users/connylin/Dropbox/Publication/Manuscript RL Alcohol hab model/Figures Tables Data/Fig1 body curve/Fig8 10sISI curve 95s/BodyCurve_10sISI/Dance_InitialEtohSensitivityPct_v1707/data.mat';
load(pData10);
D10 = MWTSet.Stats.curve.pctByPlate;
D10.ISI = repmat({'10'},size(D10,1),1);
pData60 = '/Users/connylin/Dropbox/Publication/Manuscript RL Alcohol hab model/Figures Tables Data/Fig1 body curve/Fig8 60sISI curve 95s/BodyCurve_60sISI/Dance_InitialEtohSensitivityPct_v1707/data.mat';
load(pData60);
D60 = MWTSet.Stats.curve.pctByPlate;
D60.ISI = repmat({'60'},size(D60,1),1);
% combine data
Data = [D10;D60];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170726


%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170726
data = Data.pct;
group = Data.ISI;
Nname = 'plate';
ps = fullfile(pSave,'AMOVA.txt');
[txt,anovastats,multstats,T,ST,Stats] = anova1_std_v2(data,group, Nname,'pSave',ps);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170726

%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170726
save(fullfile(pSave,'data.mat'), 'Data','Stats');
fprintf(' *** finish ***\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170726


return





