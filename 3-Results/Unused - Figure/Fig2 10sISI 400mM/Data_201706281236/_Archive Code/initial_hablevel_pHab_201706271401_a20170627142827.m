%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GET MWT INFORMATION
% get mwt paths
pDatabase = '/Users/connylin/Dropbox/Publication/Manuscript RL Alcohol hab model/Figures Tables Data/Fig1 10sISI 400mM large scale';
load(fullfile(pDatabase,'MWTDB.mat'))
pMWT = MWTDB.mwtpath;
clear MWTDB;

%% ANALYZE DATA
MWTSet = Dance_ShaneSpark4(pMWT,pM);

%% GENERATE STATS

















