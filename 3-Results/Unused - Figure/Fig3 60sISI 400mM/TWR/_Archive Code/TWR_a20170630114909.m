%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GET MWT INFORMATION
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

%% ANALYZE DATA
MWTSet = Dance_ShaneSpark_v5r1(pMWT,'pSave',pSave);




return
