%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; 
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
p = '/Users/connylin/Dropbox/RL Pub PhD Dissertation/Chapters/2-Wildtype/3-Results/Figures & Data/Fig2-1 Body Curve 60m 3 4 5 do/AssayAge/Stats_v3_selectTimes10min';
pd = fullfile(p,'curve rmANOVA.txt');

A = rmANOVAextractTxtOutputObj;
A.datapath =pd;
A.txt_t_identifier = 'Posthoc[(]Tukey[)]t by gname';

%%

%%
pvlim = 0.001;
pvsig = 0.05;
g1 = 'N2_3d';
g2 = 'N2_400mM_5d';

str = posthoc_tXgroupWriteUp(A,pvlim,pvsig);


fid = fopen(fullfile(p,'curve rmANOVA pairwise report.txt'),'w');
fprintf(fid,'%s',str);
fclose(fid);
%%
