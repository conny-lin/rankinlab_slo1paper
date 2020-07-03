%% analyze rType for all strans.
% 2020-02-16 12:15
%% add paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% set source and save paths
% p = mfilename('fullpath');
% [~,f] = fileparts(p);
% pSave = create_savefolder(p, f);
pSave = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Figures & Data/Fig5 slo1 cell specific ko';


%% define by group
StrainNameTarget = {'HKK796','HKK1165','NM1968','VG902','VG903'};


%% TWR
for i = 1:numel(StrainNameTarget)
    % get strain name
    strainname = StrainNameTarget{i};
    % display running
    fprintf('analyzing %s\n',strainname);
    
    % get MWT file paths (wildtype only from mutant exp + N = plate)-------
    % get path for TWR
    p = fullfile(pSave,strainname,'Dance_ShaneSpark_v5r1','Dance_ShaneSpark_v5r1.mat');
    load(p);
    pMWT = MWTSet.MWTDB.mwtpath;
    fprintf('got paths\n');
    
    
    % analysis - TWR
    pSaveA = fullfile(pSave,strainname);
    fprintf('analyzing TWR for %s\n',strainname);
    Dance_ShaneSpark_v5r2(pMWT,'pSave',pSaveA);

    fprintf('Done TWR ... %s\n',strainname);
 
end

% fprintf('TAR*****\n');
% 
% for i = 1:numel(StrainNameTarget)
%     % get strain name
%     strainname = StrainNameTarget{i};
%     % display running
%     fprintf('analyzing %s\n',strainname);
%     
%     % load paths to mwt files in database
%     p = fullfile(pSave,strainname,'TAR','Dance_rType','data.mat');
%     load(p)
%     pMWT = Data.mwtpath;
%     fprintf('got paths\n');
% 
%     % analysis TAR
%     fprintf('analyzing TAR for %s\n',strainname);
%     pSaveA = fullfile(pSave,strainname,'TAR');
%     Dance_rType2(pMWT,'pSave',pSaveA);
% 
%     fprintf('Done TAR %s\n',strainname);
% 
% end










