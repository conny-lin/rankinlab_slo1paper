%% analyze probabilty of acceleration and reversals
%% paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% Define save paths
pSave = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Fig2-Rev For Prob wt';

%% get experiment information 
pData = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Fig2-Rev For Prob wt/Fig2-11 response type/resp_type_prob_20170604/data.mat';
load(pData);