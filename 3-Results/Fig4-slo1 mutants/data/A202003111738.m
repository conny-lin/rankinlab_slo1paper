%% analyze rType for all strans.
% 2020-02-16 12:15
%% add paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% set source and save paths
pSave = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Figures & Data/Fig4 slo1 mutants';

%% separate by group
StrainNameTarget = {'JPS429', 'BZ142','NM1968', 'NM1630'};

for i = 1:numel(StrainNameTarget)
    % get strain name
    strainname = StrainNameTarget{i};
    % display running
    fprintf('analyzing %s\n',strainname);
    
    % load paths to mwt files in database
    p = fullfile(pSave,strainname,'TAR','Dance_rType','data.mat');
    load(p)
    pMWT = Data.mwtpath;
    fprintf('got paths\n');

    % analysis - TWR
    pSaveA = fullfile(pSave,strainname,'TWR');
    fprintf('analyzing TWR for %s\n',strainname);
    Dance_ShaneSpark_v5r1(pMWT,'pSave',pSaveA);

    fprintf('Done TWR ... %s\n',strainname);

end

fprintf('TAR*****\n');

for i = 1:numel(StrainNameTarget)
    % get strain name
    strainname = StrainNameTarget{i};
    % display running
    fprintf('analyzing %s\n',strainname);
    
    % load paths to mwt files in database
    p = fullfile(pSave,strainname,'TAR','Dance_rType','data.mat');
    load(p)
    pMWT = Data.mwtpath;
    fprintf('got paths\n');

    % analysis TAR
    fprintf('analyzing TAR for %s\n',strainname);
    pSaveA = fullfile(pSave,strainname,'TAR');
    Dance_rType2(pMWT,'pSave',pSaveA);

    fprintf('Done TAR %s\n',strainname);

end










