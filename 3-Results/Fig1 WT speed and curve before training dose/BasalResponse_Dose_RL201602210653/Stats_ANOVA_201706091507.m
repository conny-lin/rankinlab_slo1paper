%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'FB','genSave',true);
% addpath(pM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msrlist = {'speed','curve'};
incr = 5; % increment
tt = [5:incr:60]; % decide time frames to assay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pDataFolder = fullfile(fileparts(pM));
datafilename = 'data.mat';
pData = fullfile(pDataFolder,datafilename);
load(pData);
DataMaster = DataExpPCT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for msri = 1:numel(msrlist)
    
    msr = msrlist{msri}; % get msr name
    data = DataMaster.(msr); % get data
    groups = DataMaster.groupname; % get groups
    [txt,anovastats,multstats] = anova1_std_v2(data,groups,'worms'); % anova
    fid = fopen(fullfile(pM,sprintf('%s anova.txt',msr)),'w'); % open output text file name
    fprintf(fid,'%s\n',txt); % add to printout
    fclose(fid); % close file
    
end
   







