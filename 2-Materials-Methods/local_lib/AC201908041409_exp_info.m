%% analyze body curve for all strans.
% 2019-07-29-16:49
%% paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% get experiment information 
pData = '/Volumes/COBOLT/MWT';
D = load([pData,'/MWTDB.mat']);
D = D.MWTDB.text;
% get 10sISI
i = ismember(D.rc,'100s30x10s10s');
D = D(i,:);
%% get experiment with the targeted strain
strainT = {'BZ142','NM1630','NM1968','HKK1165','HKK796','VG902','VG903'};
i = ismember(D.strain,strainT);
ExpT = unique(D.expname(i));
i = ismember(D.expname,ExpT);
D = D(i,:);
% get rid of strains in exp thats not targeted strains
groupT = {'BZ142','NM1630','NM1968','HKK1165','HKK796','VG902','VG903',...
    'BZ142_400mM','NM1630_400mM','NM1968_400mM','HKK1165_400mM',...
    'HKK796_400mM','VG902_400mM','VG903_400mM','N2','N2_400mM'};
i = ismember(D.groupname,groupT);
D = D(i,:);

%% generate experiment sheet
pSave = '/Users/connylin/Dropbox/CA Publications/Manuscript RL Alcohol hab model slo1/2-Materials Methods/MWTDB slo-1.csv';

writetable(D,pSave)
