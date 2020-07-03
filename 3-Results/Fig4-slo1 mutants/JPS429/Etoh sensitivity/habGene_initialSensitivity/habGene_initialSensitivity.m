%  habituation curve 10sISI of valid experiments

%% INITIALIZING
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
pSaveM = fileparts(pM);


%% Settings
msrlist = {'RevFreq','RevSpeed','RevDur'};
pDataHome = '/Users/connylin/Dropbox/RL/PubInPrep/PhD Dissertation/4-STH genes/Data/10sIS by strains';
% Dance = 'Dance_ShaneSpark4_Nlim8';

% limit strain
strainlimit = {'BZ416'};
strainlimit = {};
analysis_done = true;
%% get strain names
strainlist = dircontent(pDataHome);
if ~isempty(strainlimit); 
    strainlist(~ismember(strainlist,strainlimit)) = []; 
end
if isempty(strainlist); warning('no strain'); end


%% run data for all strains
if ~analysis_done
for si = 1:numel(strainlist)
    % load data
    strain = strainlist{si};
    pSave = sprintf('%s/%s',pDataHome,strain);
    fprintf('%d/%d: %s **********\n',si,numel(strainlist),strain);
    clear MWTDB; 
    pD = [pSave,'/DataTrv.mat'];
    nodata = false;
    if exist(pD,'file'); load(pD,'MWTDB');
    else
       % MWTDB - get all paths
        load('/Volumes/COBOLT/MWT/MWTDB.mat','MWTDB');
        MWTDB = MWTDB.text; 
        i = MWTDB.ISI == 10 & MWTDB.preplate == 100 ...
            & MWTDB.tapN == 30 ...
            & ismember(MWTDB.groupname,sprintf('%s_400mM',strain));
        exp = unique(MWTDB.expname(i));
        MWTDB = MWTDB(ismember(MWTDB.expname,exp),:);
        MWTDB(~ismember(MWTDB.groupname,{'N2','N2_400mM',strain,[strain,'_400mM']}),:) = [];
        pMWT = MWTDB.mwtpath;
        nodata = isempty(pMWT);
    end
    if nodata; break; end
    pMWT = MWTDB.mwtpath; % get pMWT path
    mwtID = MWTDB.mwtid;
    
    
    %% settings
    pSave = create_savefolder(pSave,'InitialEtohSensitivity');
    MWTSet = struct;
    MWTSet.MWTDB = parseMWTinfo(pMWT);
    startTime = 90; 
    endTime = 95;

    %% get data
    DataMeta = IS_getData(pMWT,MWTDB.mwtid, startTime, endTime);
    
    %% stats
    IS_stats(DataMeta,MWTDB,pSave);
    
    %% calculate initial sensitivity eStats
    IS_graph(DataMeta,MWTDB,pSave,startTime,endTime);
    
    
end

end

%% collect graphs
filenames = {'Speed t90-95.pdf' 'curve MANOVA.txt' 'curve t90-95.pdf',...
    'speed MANOVA.txt' 'speedbm MANOVA.txt' 'speedbm t90-95.pdf'}';
pMF = repmat({pM}, size(filenames));


for si = 1:numel(strainlist)
    % load data
    strain = strainlist{si};
    pSave = sprintf('%s/%s',pDataHome,strain);
    fprintf('%d/%d: %s **********\n',si,numel(strainlist),strain);

    pSave = create_savefolder(pSave,'InitialEtohSensitivity');
    
    a = repmat({strain},numel(filenames),1);
    newnames = strjoinrows([a filenames]);
    
    pd = strjoinrows([pMF,newnames],'/');
    
    ps = strjoinrows([cellfunexpr(filenames,pSave) filenames],'/');
    
    cellfun(@copyfile,ps,pd)
    
end














