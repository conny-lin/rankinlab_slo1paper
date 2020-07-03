%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% GLOBAL INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths ++++++++++++++
load(fullfile(fileparts(pSave),'MWTDB.mat'),'MWTDB');
pMWT = MWTDB.mwtpath;
% ---------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIAL SENSITIVITY 
fprintf('\n\n --- CURVE ANALYSIS ---\n');
Data_Sen = Dance_InitialEtohSensitivityPct_v1707(pMWT,pSave); % run initial sensitivity
fprintf('\n\n --- CURVE ANALYSIS DONE ---\n');

return





