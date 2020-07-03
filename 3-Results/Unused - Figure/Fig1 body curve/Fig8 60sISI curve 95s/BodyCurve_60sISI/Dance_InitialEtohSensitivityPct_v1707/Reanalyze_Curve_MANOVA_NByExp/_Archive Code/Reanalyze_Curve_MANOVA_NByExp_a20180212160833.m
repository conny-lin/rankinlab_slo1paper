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


%% GROUP DATA BY EXP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D.exp_gname = strjoinrows([D.expname, D.groupname],'@');
DataByExp = statsBasicG(D.curve, D.exp_gname);
a = DataByExp.gname;
b = regexpcellout(a,'@','split');
DataByExp.expname = b(:,1);
DataByExp.gname = b(:,2);
writetable(DataByExp,fullfile(pM,'curve_NbyExp.csv'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [text,T,p,s,t,ST] = anova1_autoresults(DataByExp.mean,DataByExp.gname)
[txt,anovastats,multstats,T,ST,ST2] = anova1_std_v2(DataByExp.mean,DataByExp.gname);
fid = fopen(fullfile(pM,'ANOVA.txt'),'w');
fprintf(fid,txt);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
