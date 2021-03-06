%% analyze rType for all strans.
% 2020-02-16 12:15
%% add paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% set source and save paths
pPaperHome = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1';
pSave = fullfile(pPaperHome,'3-Results','Data_4ExpSlo1');

%% get pMWT information 
pMat = fullfile(pSave,'Dance_ShaneSpark4','Dance_ShaneSpark4.mat');
load(pMat);
pMWT = MWTSet.MWTDB.mwtpath;

%% separate by group
StrainNameTarget = {'HKK1165' 'HKK796', 'NM1968', 'VG902', 'VG903'};
[GroupedpMWT4Dance] = groupMWTpaths4Dance(pMWT,StrainNameTarget);

for i = 1:numel(StrainNameTarget)
    pMWT = GroupedpMWT4Dance{i};
    pSaveA = fullfile(pSave,'Dance_rType2');
    pSaveA1 = create_savefolder(pSaveA,StrainNameTarget{i});
    
    %% analysis
    % Dance_InitialEtohSensitivityPct(pMWT,'pSave',pSave);
    % Dance_Raster2(pMWT,'pSave',pSave);
    % 
    % Dance_ShaneSpark4(pMWT, pSave);
    % Dance_ShaneSpark_v5r1(pMWT,'pSave',pSave);
     Dance_rType2(pMWT,'pSave',pSaveA1);

end









