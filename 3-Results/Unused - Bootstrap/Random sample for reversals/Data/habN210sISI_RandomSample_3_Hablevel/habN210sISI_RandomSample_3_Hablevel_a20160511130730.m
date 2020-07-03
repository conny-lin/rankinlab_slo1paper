%  habituation curve 10sISI of valid experiments

%% INITIALIZING
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
pSave = fileparts(pM);

%% Settings
plateN = 9;
samplingN = 100;
sampleRepeatN = 3;
msrlist = {'RevFreq','RevSpeed','RevDur'};
phenotype = 'HabLevel';
postfix = sprintf('%s plateN%d sample%d',phenotype,plateN,samplingN);

%% LOAD MWT DATABASE
load([pSave,'/DataTrv.mat']);

%% get habituated taps
DataTrv(DataTrv.tap<=28,:) = [];

%% get plate mean
Data = table;
Data.mwtid = unique(DataTrv.mwtid);
for msri = 1:numel(msrlist)
    d = DataTrv.(msrlist{msri});
    T = grpstatsTable(d, DataTrv.mwtid);
    Data.(msrlist{msri}) = T.mean;
    
end


%% put in groups
[i,j] = ismember(DataTrv.mwtid,MWTDB.mwtid);
DataTrv.groupname = MWTDB.groupname(j(i));
groupNames = unique(DataTrv.groupname);
Data = cell(size(groupNames));
for gi = 1:numel(groupNames)
    Data{gi} = DataTrv(ismember(DataTrv.groupname,groupNames(gi)),:);     
end

%% random sample 
report = nan(sampleRepeatN,numel(msrlist));
for sri = 1:sampleRepeatN
    
for msri = 1:numel(msrlist)
    Mean = nan(samplingN,numel(groupNames));
    pValue = nan(samplingN,1);
    N = nan(size(Mean));
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
        % get rid of nan
        i = isnan(Y) | isinf(Y);
        if sum(i) > 0
           Y(i) = [];
           G(i) = [];
        end
        
        
        [~,T,p] = anova1_autoresults(Y,G);
        pValue(samplei) = p;
        
        [i,j] = ismember(groupNames,T.gnames);
        
        Mean(samplei,:) = T.mean(j(i));
        N(samplei,:) = T.N(j(i));

    end
    
    % export result
    MeanDiff = Mean(:,2) - Mean(:,1);
    
    T = table;
    T.(groupNames{1}) = Mean(:,1);
    T.(groupNames{2}) = Mean(:,2);
    T.Diff = MeanDiff;
    T.pValue = pValue;
    savename = sprintf('%s %s r%d.csv',msr,postfix,sri);
    cd(pSave); writetable(T,savename);
    
    % calculate percent
    nP = sum(pValue < 0.05)./samplingN;
    fprintf('%s %s: %f.0%% significant\n',msr,phenotype,nP*100)
    report(sri,msri) = nP;
end

end

%% export report

a = array2table(report,'VariableNames',msrlist);
cd(pSave); writetable(a,sprintf('%s percent significant a95p.csv',phenotype));

return






















