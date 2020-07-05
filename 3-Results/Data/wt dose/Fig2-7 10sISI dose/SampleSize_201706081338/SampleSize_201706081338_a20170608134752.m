%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; 
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
load(fullfile(fileparts(pM),'MWTDB.mat'))

A = tabulate(MWTDB.groupname);
n = cell2mat(A(:,2));
str = gen_Nstring(n);

fid = fopen(fullfile(pM,'N.txt'));
fprintf('%s',str)
fclose(fid);