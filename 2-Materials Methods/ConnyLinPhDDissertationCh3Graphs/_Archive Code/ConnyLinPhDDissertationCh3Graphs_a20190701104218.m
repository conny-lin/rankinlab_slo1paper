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
%% PROGRAMS (DONE)
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MWTDatabase_v2(pNewData,pDataBase);
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% program (current)
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% load MWTDB
p = sprintf('%s/%s',pDataBase,'MWTDB.mat');
load(p);
clear p;

%% find slo-1 strains -----------------------------------------------------
a = MWTDB.strain.genotype;
searchTerm = 'slo-1';
c = false(size(a));
for i = 1:numel(a)
    b = regexp(a{i},searchTerm,'once');
    if ~isempty(b)
        c(i) = true;
    end
    clear b;
end
% report result
n = sum(c);
fprintf('# of strains found: %d\n',n);
% disp(a(c))
% get table
StrainInfo = MWTDB.strain(c,:);
% save slo-1 strain information in csv
p = sprintf('%s/2-Materials Methods/Strain Info - slo-1.csv',pPaper);
writetable(StrainInfo,p);
clear a i searchTerm c n p;

%% get how many plates per strain
DB = MWTDB.text; % get experiment database
strains = StrainInfo.strain; % get strain names
[i,j] = ismember(StrainInfo.strain, DB.strain); % find mwt plates with slo-1 strains
iSLO1 = i;
% find experiment name unique
expnameU = unique(DB.expname(iSLO1));
n = numel(expnameU);
fprintf('# of unique experiments found: %d\n',n);


%% find 10sISI 
rc = DB.rc;
rcu = unique(rc);
n = numel(rcu);
fprintf('# of unique runing conditions found: %d\n',n);
rcu(iSLO1)


return

%--------------------------------------------------------------------------
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% program (to do)
MWTSet = Dance_TWR1(pMWT,pSave);







































