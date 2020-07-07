%% global setting
filesep = '/'; % change for PC

%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% input paths (testing)
pNewData = '/Users/connylin/Dropbox/RL Data new 20190608/Conny_rawdata for analysis copy 2';
pDataBase = '/Users/connylin/Dropbox/RL Data new 20190608/MWT test database';
pSave = '/Users/connylin/Dropbox/RL Data new 20190608/Output';

%% input paths (Conny's default)
% pNewData = '/Volumes/COBOLT/MWT_New';
% pDataBase = '/Volumes/COBOLT/MWT';

%% paths (program)
p = mfilename('fullpath');
p1 = fileparts(p);
addpath_allsubfolders(p1);
clear functionfoldername p1 p;

%% programs (done)
MWTDB = MWTDatabase_v2(pNewData,pDataBase);

return
%% program (current)
DbT = MWTDatabase_query(pData,varargin)

%% program (to do)
MWTSet = Dance_TWR1(pMWT,pSave);

















