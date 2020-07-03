%  habituation curve 10sISI of valid experiments

%% INITIALIZING
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
pSave = fileparts(pM);

%% Settings
plateN = 9;
samplingN = 100;
msrlist = {'RevFreq','RevSpeed','RevDur'};

%% LOAD MWT DATABASE
load([pSave,'/DataTrv.mat']);


%% get first tap
DataTrv(DataTrv.tap~=1,:) = [];
%% put in groups
[i,j] = ismember(DataTrv.mwtid,MWTDB.mwtid);
DataTrv.groupname = MWTDB.groupname(j(i));
groupNames = unique(DataTrv.groupname);
Data = cell(size(groupNames));
for gi = 1:numel(groupNames)
    Data{gi} = DataTrv(ismember(DataTrv.groupname,groupNames(gi)),:);     
end

%% random sample 
for msri = 1:numel(msrlist)
    Stats
    Mean = nan(samplingN,numel(groupNames));
    pValue = nan(samplingN,1);
    
    for samplei =  1:samplingN
        Y = nan(plateN,numel(groupNames));
        G = reshape(repmat(groupNames',plateN,1),numel(Y),1);

        for gi = 1:numel(groupNames)
            gn = groupNames{gi};
            msr = msrlist{msri};
            data = Data{gi}.(msr);
            Y(:,gi) = datasample(data,plateN);
        end
        
        Y = reshape(Y, numel(Y),1);
        
        [~,T,p] = anova1_autoresults(Y,G);
        pValue(samplei) = p;
        [i,j] = ismember(groupNames,T.gnames);
        Mean(samplei,:) = T.mean(j(i));
        

    end
    
    return
end
return






















