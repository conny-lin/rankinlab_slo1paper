%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; 
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data = load(fullfile(fileparts(pM),'data.mat'));

D = Data.DataMeta;
D = innerjoin(D,Data.MWTDB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GROUP DATA BY PLATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D.exp_gname = strjoinrows([D.expname, D.groupname],'@');
DataByPlate = statsBasicG(D.curve, D.mwtid);

return
a = DataByPlate.gname;
b = regexpcellout(a,'@','split');
DataByPlate.expname = b(:,1);
DataByPlate.gname = b(:,2);
writetable(DataByPlate,fullfile(pM,'curve_NbyExp.csv'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [text,T,p,s,t,ST] = anova1_autoresults(DataByExp.mean,DataByExp.gname)
[txt,anovastats,multstats,T,ST,ST2] = anova1_std_v2(DataByPlate.mean,DataByPlate.gname);
fid = fopen(fullfile(pM,'ANOVA.txt'),'w');
fprintf(fid,txt);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
