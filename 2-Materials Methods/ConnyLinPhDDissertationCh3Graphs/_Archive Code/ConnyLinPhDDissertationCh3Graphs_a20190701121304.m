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
pSave = fileparts(pM);

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

%% query conditions
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
    };
condRC = '100s30x10s10s';

%% for each strains, find the experiment sets
DB = MWTDB.text;
% keep only 10s ISI
DB(~ismember(DB.rc,'100s30x10s10s'),:) = []; 
% get exp for each strain 
ExpSets = cell(size(condStrains));
eu = unique(DB.expname); % unique experiment list
for si = 1:numel(condStrains)
    strainname = condStrains{si};
    condGN = {'N2','N2_400mM',strainname,sprintf('%s_400mM',strainname)};
    % set qualifying matrix for unique experiment list
    m = false(size(eu,1),numel(condGN));
    for gi = 1:numel(condGN)
        % index to mwt with the group name
        i = ismember(DB.groupname,condGN(gi));
        % find experiment that has the group name
        u = unique(DB.expname(i));
        % index to experiment list
        j = ismember(eu,u);
        % mark that experiment as containing the group
        m(j,gi) = true;
    end
    % get index to qualified experiments on the list
    i = sum(m,2) == numel(condGN);
    % keep only those experiments
    j = ismember(DB.expname,eu(i));
    D = DB(j,:);
    % remove plates that do not have the targeted group names
    i = ismember(D.groupname,condGN);
    D(~i,:) = [];
    % put the dataset in the file
    ExpSets{si} = D;
    % export in csv
    filename = sprintf('%s/%s_experiment info.csv',pSave,strainname);
    writetable(D,filename);
    return
end

return

i = ismember(DB.strain,condStrains); 
DB(~i,:) = [];

% query - rc
i = ismember(DB.rc,condRC);
DB(~i,:) = [];







return

%--------------------------------------------------------------------------
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% program (to do)
MWTSet = Dance_TWR1(pMWT,pSave);







































