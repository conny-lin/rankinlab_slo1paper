% convert shankspark .mat data to csv
% add function path
addpath('/Users/connylin/Code/language/matlab_lib/general')
% define data dir
data_dir = '/Users/connylin/Dropbox/CA/_Publications/Manuscript RL Alcohol hab model slo1/rankinlab_slo1paper/3-Results/Fig5-slo1 cell specific ko/data';
% get strain folders
strain_name = dircontent(data_dir);
% loop through strains
for folderi = 1:numel(strain_name)
    % get path to mat
    strain_dir = fullfile(data_dir,strain_name{folderi},'Dance_ShaneSpark4');
    % load mat
    load(fullfile(strain_dir,'Dance_ShaneSpark4.mat'));
    % 2. convert raw to csv
    T = MWTSet.Raw;
    writetable(T, fullfile(strain_dir, 'rawdata.csv'));
    % 3. convert db to csv
    T = MWTSet.MWTDB;
    writetable(T, fullfile(strain_dir, 'mwtdb.csv'));  
end
