%% global setting
filesep = '/'; % change for PC

%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% input paths (testing)
% pNewData = '/Users/connylin/Dropbox/RL Data new 20190608/Conny_rawdata for analysis copy 2';
% pDataBase = '/Users/connylin/Dropbox/RL Data new 20190608/MWT test database';
pSave = '/Users/connylin/Dropbox/RL Data new 20190608/Output';

%% input paths (Conny's default)
pNewData = '/Volumes/COBOLT/MWT_New';
pDataBase = '/Volumes/COBOLT/MWT';
pSave = '/Users/connylin/Dropbox/RL Data new 20190608/Output';


%% paths (program)
p = mfilename('fullpath');
p1 = fileparts(p);
addpath_allsubfolders(p1);
clear functionfoldername p1 p;

%% programs (done)
% MWTDatabase_v2(pNewData,pDataBase);


%% program (current)
p = sprintf('%s/%s',pDataBase,'MWTDB.mat');
load(p);
clear p;

return
expnames = ...
    {'20190127X_XX_100s30x10s10s_slo1';
    '20190330X_XX_100s30x10s10s_slo1';
    '20190331X_XX_100s30x10s10s_slo1';
    '20190418X_XX_100s30x10s10s_slo1'
    };
conditions = {'exp_date', expnames};

DbT = MWTDatabase_query2(pDataBase,conditions);
clear conditions

%% MWTDatabase_query
% query MWTDatabase
% Input: 
%         pData = '/Volumes/COBOLT/MWT';
%     output settings
%         DisplayVar - variables asked to display on screen in cell arrays
%     Query terms
%         mwt_id   
%         mwt      
%         mwtpath  
%         expname  
%         exp_date 
%         tracker  
%         expter   
%         rc       
%         groupname = {'N2'};
%         strain   
%         rx       
%         preplate = 100;
%         ISI = 10;
%         postrec  
%         tapN = 30;
%         genotype
% 
%     Query selection terms
%         controlgroup = {'N2'};
%         gnameSearchType:
%             'any' (default) - groups can appear in any experiment
%             'within' - all groups mentioned must be within the same experiment 
%             'withcontrol' - groups must show up with controlgroup specified in controlgroup

return

%% program (to do)
MWTSet = Dance_TWR1(pMWT,pSave);

















