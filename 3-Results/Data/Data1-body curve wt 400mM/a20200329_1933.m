%% paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')
pSave = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Data1-body curve wt 400mM';

%% get experiment information 
p = fullfile(pSave,'data.mat');
load(p)
pMWT = Data.mwtpath;

%% analysis
Dance_InitialEtohSensitivityPct(pMWT,'pSave',pSave);
