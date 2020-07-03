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
    pDataH = fullfile(fileparts(pSave),strainname); % get original data path

    % get mwt paths
    A = load(fullfile(pData,strainname,'MWTDB.mat'));
    pMWT = A.MWTDB.mwtpath; clear A;
    %% TAR
    MWTSet = Dance_rType_v1707(pMWT,pSave);
    
%     
%     return
%     %% INITIAL SENSITIVITY 
%     fprintf('\n\n --- CURVE ANALYSIS ---\n');
%     Data_Sen = Dance_InitialEtohSensitivityPct_v1707(pMWT,pSave); % run initial sensitivity
%     fprintf('\n\n --- CURVE ANALYSIS DONE ---\n');
%     
%     
%     %% TWR 
%     MWTSet = Dance_ShaneSpark_v1707(pMWT,pSave);



end
return
