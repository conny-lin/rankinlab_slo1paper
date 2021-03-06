%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUTS
pData = fullfile(pSave,'Data');
strainlist = dircontent(pData);

for strainlisti = 1:numel(strainlist)
    
    strainname = strainlist{strainlisti}; % get strain name
%     pDataH = fullfile(fileparts(pSave),strainname); % get original data path
    pSaveS = fullfile(pSave,'Data',strainname); 
    % get mwt paths
    A = load(fullfile(pSaveS,'MWTDB.mat'));
    pMWT = A.MWTDB.mwtpath; clear A;
    % TWR 
    MWTSet = Dance_ShaneSpark_v1707(pMWT,pSaveS);
    
    % TAR
%     MWTSet = Dance_rType_v1707(pMWT,pSave);
    
%     
%     %% INITIAL SENSITIVITY 
%     fprintf('\n\n --- CURVE ANALYSIS ---\n');
%     Data_Sen = Dance_InitialEtohSensitivityPct_v1707(pMWT,pSave); % run initial sensitivity
%     fprintf('\n\n --- CURVE ANALYSIS DONE ---\n');
%     
%     



end
return
