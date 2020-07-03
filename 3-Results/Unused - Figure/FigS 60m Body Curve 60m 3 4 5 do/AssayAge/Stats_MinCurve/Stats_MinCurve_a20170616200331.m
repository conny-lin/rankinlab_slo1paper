%% OBJECTIVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stats for this statement:
% For worms on 400mM ethanol, 3 and 4 days old worms had a minimum curve of
% 15 degrees, while 5 days old worms had a minimum curve of 19 degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; 
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% universal settings
pData = '/Users/connylin/Dropbox/RL/RL Pub PhD Dissertation/Chapters/2-Wildtype/3-Results/1-Figures & Data/Fig2-1 Body Curve 60m 3 4 5 do/AssayAge/Dance_DrunkMoves'; % paths
varname = 'curve'; % variable name
timelist = [1, 5:5:60];

%% get data 
datafilename = 'Dance_DrunkMoves.mat'; % get data file name
fname = fullfile(pData,datafilename); % get path to file
load(fname); % load all data
% get data of interest
v = {'groupname','mwtname','timeind','mean'}; % define table variable of interest
Data = MWTSet.Data_Plate.(varname)(:,v); %% get curve data per plate with variable of interest
%% translate groupname from numeric to text
Data(~ismember(Data.timeind,timelist),:) = []; % get only time of interest
g = MWTSet.Info.VarIndex.groupname; % get groupname reference
Data.groupname = g(Data.groupname); % replace numeric index to reference

return

%% Stats
% take min curve per plate
% anova comparison 
% post hoc comparison
return
Data = readtable(fname); % import csv as table
timelist = [5:5:60]; % designate selected time
Data = Data(:,[1, timelist+1]); % get only selected time

%% find max minutes per group
D = table2array(Data(:,2:end));
rowN = size(D,1);
tmax = nan(rowN,1);
tmaxnum = tmax;
tmin = tmax;
tminnum = tmax;
for ri = 1:rowN
    [m,t] = max(D(ri,:)); % find max per row
    tmax(ri) = timelist(t); % translate to time
    tmaxnum(ri) = m;

    [m,t] = min(D(ri,:)); % find max per row
    tmin(ri) = timelist(t); % translate to time
    tminnum(ri) = m;
end

T = table;
T.groupname = Data.groupname;
T.timemax = tmax;
T.(sprintf('%s_max',varname)) = tmaxnum;
T.timemin = tmin;
T.(sprintf('%s_min',varname)) = tminnum;

filename = fullfile(pM,sprintf('%s max.csv',varname));
writetable(T,filename)

