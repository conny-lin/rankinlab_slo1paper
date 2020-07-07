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
c = false(size(a));
for i = 1:numel(a)
    b = regexp(a{i},'slo-1','once');
    if ~isempty(b)
        c(i) = true;
    end
    clear b;
end
% report result
n = sum(c);
fprintf('# of strains found: %d\n',n);
disp(a(c))
% get table
StrainInfo = MWTDB.strain(c,:);
% save slo-1 strain information in csv
p = sprintf('%s/2-Materials Methods/Strain Info - slo-1.csv',pPaper);
writetable(StrainInfo,p);
clear a i searchTerm c n p;

%% get how many plates per strain
DB = MWTDB.text; % get experiment database
strains = StrainInfo.strain; % get strain names
i = ismember(DB.strain,StrainInfo.strain); % find mwt plates with slo-1 strains
DB(~i,:) = []; % delete ones without slo-1
% get 10sISI, slo-1 exp
i = ismember(DB.rc,'100s30x10s10s');
DB(~i,:) = [];
% get unique exp name
enu = unique(DB.expname);

%% get data base only with the same exp name
DB = MWTDB.text;
i = ismember(DB.expname,enu);
DB(~i,:) = [];
size(DB)



return

%--------------------------------------------------------------------------
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% program (to do)
MWTSet = Dance_TWR1(pMWT,pSave);







































