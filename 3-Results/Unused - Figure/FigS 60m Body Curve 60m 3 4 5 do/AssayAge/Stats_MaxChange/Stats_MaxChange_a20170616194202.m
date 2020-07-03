%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; 
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get paths
pData = '/Users/connylin/Dropbox/RL/RL Pub PhD Dissertation/Chapters/2-Wildtype/3-Results/1-Figures & Data/Fig2-1 Body Curve 60m 3 4 5 do/AssayAge/Dance_DrunkMoves/Data Curve';

%% load data
varnames_list = {'curve','speed'};
for vi = 1:numel(varnames_list)
    varname = varnames_list{vi};
    fname = fullfile(pData,sprintf('%s.csv',varname)); % file name
    Data = readtable(fname); % import csv as table
    timelist = [5:5:60]; % designate selected time
    Data = Data(:,[1, timelist+1]); % get only selected time

    %% find max minutes per group
    D = table2array(Data(:,2:end));
    rowN = size(D,1);
    tmax = nan(rowN,1);
    tmaxnum = tmax;
    for ri = 1:rowN
        [m,t] = max(D(ri,:)); % find max per row
        tmax(ri) = timelist(t); % translate to time
        tmaxnum(ri) = m;
    end

    T = table;
    T.groupname = Data.groupname;
    T.timemax = tmax;
    T.(varname) = tmaxnum;

    filename = fullfile(pM,sprintf('%s max.csv',varname));
    writetable(T,filename)

end