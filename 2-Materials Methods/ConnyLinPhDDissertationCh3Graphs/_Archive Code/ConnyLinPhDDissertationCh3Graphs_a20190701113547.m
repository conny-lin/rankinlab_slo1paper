%% Information

%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% input paths (Conny's default)
pNewData = '/Volumes/COBOLT/MWT_New';
pDataBase = '/Volumes/COBOLT/MWT';
pPaper = '/Users/connylin/Dropbox/CA Publications/Manuscript RL Alcohol hab model slo1';
pSave = '/Users/connylin/Dropbox/RL Data new 20190608/Output';

%% paths (program)
p = mfilename('fullpath');
p1 = fileparts(p);
addpath_allsubfolders(p1);
clear functionfoldername p1 p;

%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% QUERY DATABASE
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% conditions
condStrains = {...
    'BZ142'
    'N2'
    'NM1630'
    'NM1968'
    'JPS428'
    'JPS429'
    'HKK1165'
    'HKK796'
    'VG902'
    'VG903'        
    'BZ416'
    'CX3940'
    };
condRC = '100s30x10s10s';
    
% get database
p = sprintf('%s/%s',pDataBase,'MWTDB.mat');
load(p);
clear p;
DB = MWTDB.text; % get experiment database

% query - strain
i = ismember(DB.strain,condStrain); 
DB(~i,:) = []; 
% query - rc
i = ismember(DB.rc,condRC);
DB(~i,:) = [];




return

%--------------------------------------------------------------------------
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% program (to do)
MWTSet = Dance_TWR1(pMWT,pSave);







































