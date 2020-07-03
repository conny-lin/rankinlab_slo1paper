%% OBJECTIVE:
% - create a list of experiments relevant for rapid tolerance chapter
%% INITIALIZING
clc; clear; close all;
%% PATHS
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);

%% load database
pDB = '/Volumes/COBOLT/MWT/MWTDB.mat';
load(pDB);
MWTDB = MWTDB.text;

%% target list
rclist = {'3600s0x0s0s'};

%% query
i = ismember(MWTDB.rc, rclist) & ismember(MWTDB.strain,'N2')...
    & ~ismember(MWTDB.groupname, {'N2_Test','N2_400mM_NoFood','N2_NoFood'});
MWTDB(~i,:) = [];
%% change to group name alternate
cd('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/MWTDatabase')
rxname = readtable('rxname_alternate.csv');
for ri =1:size(rxname,1)
    i= ismember(MWTDB.rx,rxname.rx(ri));
    MWTDB.rx(i) = rxname.rx_alternate(ri);
end
% create new names
A = [MWTDB.strain MWTDB.rx];
MWTDB.groupname = strjoinrows(A,'_');

%% keep only 1 hour recovery time
% Rx = parseRx(MWTDB.rx);
% MWTDB(Rx.rec_hr~=1,:) = [];
unique(MWTDB.groupname)


%%
pMWT = MWTDB.mwtpath;
MWTSet = Dance_DrunkMoves(pMWT,'pSave',pSave);













































