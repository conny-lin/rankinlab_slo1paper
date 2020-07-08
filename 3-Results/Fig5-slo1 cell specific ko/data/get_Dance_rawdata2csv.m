%% get Dance rtype2 rawdata
%  Conny Lin | June 8, 2020
%% add paths (program)
addpath(genpath('/Users/connylin/Code/language/matlab_lib'));
addpath(genpath('/Users/connylin/Code/proj/rankin_lab'));
% set source and save paths
pSave = '/Users/connylin/Dropbox/CA/_Publications/Manuscript RL Alcohol hab model slo1/rankinlab_slo1paper/3-Results/Fig5-slo1 cell specific ko/data';
[~,~,strain_names, strain_folders] = dircontent(pSave);

%% getdata / rType2
for i = 1:numel(strain_names)
    % get strain specific variables
    strainname = strain_names{i};
    pdata = fullfile(pSave, strainname, 'Dance_rType2', 'data.mat');
    % load data
    load(pdata)
    % save csv
    writetable(Data, ...
        fullfile(pSave, strainname, 'Dance_rType2','rawdata.csv'))
end

%% get data / ShaneSpark_v5r1
for i = 1:numel(strain_names)
    % get strain specific variables
    strainname = strain_names{i};
    pdata = fullfile(pSave, strainname, 'Dance_ShaneSpark_v5r1', 'Dance_ShaneSpark_v5r1.mat');
    % load data
    load(pdata)
    % save mwtdb
    writetable(MWTSet.MWTDB, ...
        fullfile(pSave, strainname, 'Dance_ShaneSpark_v5r1','mwtdb.csv'));
    % save raw data
    writetable(MWTSet.Raw, ...
        fullfile(pSave, strainname, 'Dance_ShaneSpark_v5r1','rawdata.csv'));
end


