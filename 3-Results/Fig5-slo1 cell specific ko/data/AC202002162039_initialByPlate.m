%% get data by plate for initial response
% 2020-02-16 20:39
%% add paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% set source and save paths
pHome = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Data_4ExpSlo1/AC20200216';


%% get strain names
[~,~,fn,~] = dircontent(pHome);

%% go through each strain
for fni = 1:numel(fn)
    
    pData = fullfile(pHome,fn{fni},'Dance_InitialEtohSensitivityPct','data.mat');
    DataByPlate = initialResponseByPlate(pData);
    
end