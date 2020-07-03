%  habituation curve 10sISI of valid experiments

%% INITIALIZING
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
pSave = fileparts(pM);

%% LOAD MWT DATABASE
load([pSave,'/MWTInfo.mat']);
pMWT = MWTDB.mwtpath;
%% rename MWTDB.mwtid
MWTDB.mwtid = [1:size(MWTDB,1)]';
%% INITIAL RESPONSE
DataTrv = import_trv_table(pMWT);
%%
cd(pSave); save('DataTrv.mat','DataTrv','MWTDB');

return










