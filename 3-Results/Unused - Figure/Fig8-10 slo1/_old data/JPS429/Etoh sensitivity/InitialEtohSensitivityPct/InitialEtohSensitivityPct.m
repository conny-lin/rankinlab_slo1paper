%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% GLOBAL INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths ++++++++++++++
pData_Strain = '/Users/connylin/Dropbox/RL Pub PhD Dissertation/Chapters/3-Genes/3-Results/0-Data/10sIS by strains';
% --------------------

% settings ++++++
pvsig = 0.05;
pvlim = 0.001;
startTime = 90; 
endTime = 95;
msrlist = {'RevFreq','RevSpeed','RevDur'};
% strains
strainNames = DanceM_load_strainInfo;
% ---------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings


%% run data for all strains
for si = 1:size(strainNames,1)
    
    % get strain +++
    processIntervalReporter(size(strainNames,1),1,'*** strain',si);
    strain = strainNames.strain{si};
    genotype = strainNames.genotype{si};
    % --------------
    
    % get save path +++
    pSave = create_savefolder(fullfile(pData_Strain,strain,'Etoh sensitivity',mfilename));
    % copy code to folder
    source = [pM,'.m'];
    destination = fullfile(pSave,[mfilename,'.m']);
    copyfile(source,destination)
    % -----------------

    % load data ++++++++
    pD = fullfile(pData_Strain,strain,'TWR/DataTrv.mat');
    load(pD,'MWTDB');
%    % MWTDB - get all paths
%     load('/Volumes/COBOLT/MWT/MWTDB.mat','MWTDB');
%     MWTDB = MWTDB.text; 
%     i = MWTDB.ISI == 10 & MWTDB.preplate == 100 ...
%         & MWTDB.tapN == 30 ...
%         & ismember(MWTDB.groupname,sprintf('%s_400mM',strain));
%     exp = unique(MWTDB.expname(i));
%     MWTDB = MWTDB(ismember(MWTDB.expname,exp),:);
%     MWTDB(~ismember(MWTDB.groupname,{'N2','N2_400mM',strain,[strain,'_400mM']}),:) = [];
%     pMWT = MWTDB.mwtpath;
%     nodata = isempty(pMWT);
% 
%     pMWT = MWTDB.mwtpath; % get pMWT path
%     mwtID = MWTDB.mwtid;
%     
    % -------------------

    
%     MWTSet = struct;
%     MWTSet.MWTDB = parseMWTinfo(pMWT);
    %% get data from gangnam ++++++
    DataMeta = IS_getData(MWTDB.mwtpath,MWTDB.mwtid, startTime, endTime);
    % ----------------------------
    
    %% calculate percentage difference +++++++++
    
    
    % -----------------------------------------
    
    
    
    %% stats
    IS_stats(DataMeta,MWTDB,pSave);
    
    %% calculate initial sensitivity eStats
    IS_graph(DataMeta,MWTDB,pSave,startTime,endTime);
    
    % save +++
    save(fullfile(pSave,'data.mat'),'DataMeta','MWTDB')
    % -------
end


return
%% collect graphs
% filenames = {'Speed t90-95.pdf' 'curve MANOVA.txt' 'curve t90-95.pdf',...
%     'speed MANOVA.txt' 'speedbm MANOVA.txt' 'speedbm t90-95.pdf'}';
% pMF = repmat({pM}, size(filenames));
% 
% 
% for si = 1:size(strainNames,1)
%     % load data
%     strain = strainNames.strain{si};
%     pSave = sprintf('%s/%s',pData_Strain,strain);
%     fprintf('%d/%d: %s **********\n',si,numel(strainlist),strain);
% 
%     pSave = create_savefolder(pSave,'InitialEtohSensitivity');
%     
%     a = repmat({strain},numel(filenames),1);
%     newnames = strjoinrows([a filenames]);
%     
%     pd = strjoinrows([pMF,newnames],'/');
%     
%     ps = strjoinrows([cellfunexpr(filenames,pSave) filenames],'/');
%     
%     cellfun(@copyfile,ps,pd)
%     
% end














