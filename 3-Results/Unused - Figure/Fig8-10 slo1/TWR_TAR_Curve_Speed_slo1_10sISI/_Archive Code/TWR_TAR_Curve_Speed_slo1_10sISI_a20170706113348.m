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
    pSaveS = fullfile(pSave,'Data',strainname); % save path
    A = load(fullfile(pSaveS,'MWTDB.mat')); % get mwt paths
    pMWT = A.MWTDB.mwtpath; clear A; % pMWT
    
%     MWTSet = Dance_ShaneSpark_v1707(pMWT,pSaveS); % TWR 
    
    MWTSet = Dance_rType_v1707(pMWT,pSave); % TAR
    
%     
%     Data_Sen = Dance_InitialEtohSensitivityPct_v1707(pMWT,pSaveS); % run initial sensitivity

return

end

