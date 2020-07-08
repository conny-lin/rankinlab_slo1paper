%% run standard analysis suites
%  Conny Lin | June 8, 2020
%% add paths (program)
addpath(genpath('/Users/connylin/Code/language/matlab_lib'));
addpath(genpath('/Users/connylin/Code/proj/rankin_lab'));

%% set source and save paths
pSave = '/Users/connylin/Dropbox/CA/_Publications/Manuscript RL Alcohol hab model slo1/rankinlab_slo1paper/3-Results/Fig5-slo1 cell specific ko/data';
[~,~,strain_names, strain_folders] = dircontent(pSave);

%% analysis
for i = 1:numel(strain_names)
    % get strain specific variables
    strainname = strain_names{i};
    psavestrain = fullfile(pSave, strainname);
    % display running
    fprintf('analyzing %s\n', strainname);
    % mwtpath
    pMWT = readtable(fullfile(strain_folders{i}, 'mwtpath.csv'), 'Delimiter',',');
    pMWT = table2cell(pMWT);
    % initial
    Dance_InitialEtohSensitivityPct(pMWT,'pSave',psavestrain);
%     pData_initial = fullfile(pSaveA, 'Dance_InitialEtohSensitivityPct','data.mat');
%     DataByPlate = initialResponseByPlate(pData_initial);
    % TWR
    Dance_ShaneSpark4(pMWT, psavestrain);
    Dance_ShaneSpark_v5r1(pMWT,'pSave',psavestrain, 'ctrl', 'HKK796');
    % TAR
    Dance_rType2(pMWT,'pSave', psavestrain);
end










