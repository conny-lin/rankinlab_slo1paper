%% analyze body curve for all strans.
% 2019-07-29-16:49
%% paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')
pSave = '/Users/connylin/Dropbox/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Data_4ExpSlo1';

%% get experiment information 
pData = '/Volumes/COBOLT/MWT';
D = load([pData,'/MWTDB.mat']);
D = D.MWTDB.text;

% get experiment with the targeted strain
expT = {'20190127X_XX_100s30x10s10s_slo1'; 
'20190330X_XX_100s30x10s10s_slo1'; 
'20190331X_XX_100s30x10s10s_slo1'; 
'20190418X_XX_100s30x10s10s_slo1'};
i = ismember(D.expname,expT);
D = D(i,:);


%%
writetable(D,fullfile(pSave,'MWTDB.csv'));