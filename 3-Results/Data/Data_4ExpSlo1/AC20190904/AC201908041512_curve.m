%% analyze body curve for all strans.
% 2019-07-29-16:49
%% paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')
pSave = '/Users/connylin/Dropbox/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Data_4ExpSlo1';

%% get experiment information 
D = readtable(fullfile(pSave,'MWTDB.csv'));
pMWT = D.mwtpath;

%% analysis
% Dance_InitialEtohSensitivityPct(pMWT,'pSave',pSave);
% Dance_Raster2(pMWT,'pSave',pSave);

Dance_ShaneSpark4(pMWT, pSave);
Dance_ShaneSpark_v5r1(pMWT,'pSave',pSave);
Dance_rType2(pMWT,'pSave',pSave);
