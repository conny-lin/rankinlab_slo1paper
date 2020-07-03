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
% get database
p = sprintf('%s/%s',pDataBase,'MWTDB.mat');
load(p);
clear p;
DB = MWTDB.text; % get experiment database

% query - strain
condStrains = {...
    'BZ142'
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
i = ismember(DB.strain,condStrains); 
DB(~i,:) = []; 

% query - rc
condRC = '100s30x10s10s';
i = ismember(DB.rc,condRC);
DB(~i,:) = [];

%% query - must have N2 and N2_400mM controls in the experiment
N2Ctrl = {'N2','N2_400mM'};
m = false(size(DB,1),numel(N2Ctrl));
for a = 1:numel(N2Ctrl)
    i = ismember(DB.groupname,N2Ctrl);
    m(i,a) = true;
end
j = sum(m,2)==numel(N2Ctrl); % must satisfy all criteria
enu = unique(DB.expname(j));
i = ismember(DB.expname,enu);
sum(i)





return

%--------------------------------------------------------------------------
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% program (to do)
MWTSet = Dance_TWR1(pMWT,pSave);







































