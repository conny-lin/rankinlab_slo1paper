%% analyze body curve for all strans.
% 2019-07-29-16:49
%% paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% testing paths
pNewData = '/Users/connylin/Dropbox/RL Data new 20190608/Conny_rawdata for analysis copy 2';
pDataBase = '/Users/connylin/Dropbox/RL Data new 20190608/MWT test database';
pSave = '/Users/connylin/Dropbox/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Data_4ExpSlo1';
%% get exp info
strainT = {'BZ142','NM1630','NM1968','HKK1165','HKK796','VG902','VG903'};

pDataBase = fullfile(pSave,'MWTDB slo-1.csv');
T = readtable(pDataBase);
pMWT = T.mwtpath;

%% body curve
[~,pStrain] = dircontent(pSH);

for si = 1:numel(pStrain)    
    pD = fullfile(pStrain{si},'Dance_rType2','data.mat');
    load(pD)
    pMWT = MWTSet.PATHS.pMWT;   
    MWTSet = Dance_InitialEtohSensitivityPct(pMWT,'pSave',pStrain{si});
    
end




