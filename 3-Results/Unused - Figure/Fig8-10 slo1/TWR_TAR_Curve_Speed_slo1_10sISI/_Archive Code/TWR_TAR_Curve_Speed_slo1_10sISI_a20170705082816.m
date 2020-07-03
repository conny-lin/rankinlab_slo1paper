%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GET MWT INFORMATION

strains = dircontent(fullfile(fileparts(pSave),'Data'));
for si = 1:numel(strains)
    pd = fullfile(pSave,'Data',strains{si});
    pd = create_savefolder(pd);
    ps = fullfile(pSave,'Data',strains{si},'MWTDB.mat');
    
    copyfile(ps,pd)
    
    return
end


return
% get mwt paths
pDatabase = fullfile('/Volumes/COBOLT/MWT','MWTDB.mat');
load(pDatabase);
MWTDB = MWTDB.text;
% select data
M = MWTDB(...
        ismember(MWTDB.rc,'100s30x60s10s') & ...
        ismember(MWTDB.groupname,{'N2','N2_400mM'})...
        ,:);
pMWT = M.mwtpath;
clear MWTDB;

%% TWR 
%% ANALYZE DATA
MWTSet = Dance_ShaneSpark_v5r1(pMWT,'pSave',pSave);


%% INITIAL SENSITIVITY 
fprintf('\n\n --- CURVE ANALYSIS ---\n');
Data_Sen = Dance_InitialEtohSensitivityPct_v1707(pMWT,pSave); % run initial sensitivity
fprintf('\n\n --- CURVE ANALYSIS DONE ---\n');



return
