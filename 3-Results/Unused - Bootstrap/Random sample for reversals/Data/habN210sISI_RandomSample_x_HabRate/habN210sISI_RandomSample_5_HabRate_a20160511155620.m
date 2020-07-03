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
phenotype = 'HabPct';
postfix = sprintf('%s plateN%d sample%d',phenotype,plateN,samplingN);

%% LOAD MWT DATABASE
load([pSave,'/DataTrv.mat']);

%% get habituation level
DataInput = DataTrv(DataTrv.tap >=28,:);
Data = table;
Data.mwtid = unique(DataInput.mwtid);
for msri = 1:numel(msrlist)
    d = DataInput.(msrlist{msri});
    T = grpstatsTable(d, DataInput.mwtid);
    Data.(msrlist{msri}) = T.mean;
    
end
HabLevel = Data;

%% get initial
Initial = DataTrv(DataTrv.tap ==1,:);

%% get asymtotic level (within experiments)
models = {'power1','power2','exp1','smoothingspline'};
mwtiu = unique(DataTrv.mwtid);
colNames = {'mwtid','f','gof_r','gof_sse'};

% Set up fittype and options.
ft = fittype( 'rat11' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [0.15 4.5 0.7];

for msri = 1:numel(msrlist)
    
   fitObj = struct;
   fitObj.mwtid = mwtiu;
   fitObj.f = cell(numel(mwtiu),1);
   fitObj.rsquare = nan(numel(mwtiu),1);
   fitObj.sse = nan(numel(mwtiu),1);
   for mwti = 3%:numel(mwtiu)
        i = DataTrv.mwtid==mwtiu(mwti);
        D = DataTrv(i,:);
        msr = msrlist{msri};
        % get data
        d = D.(msr);
        t = D.tap;

        % Fit model to data.
        [f, gof,o] = fit(t,d,ft,opts);
        plot(f,t,d)
        
        % store objects
        fitObj.mwtid(mwti) = mwtiu(mwti);
        fitObj.f{mwti} = f;
        fitObj.rsquare(mwti) = gof.rsquare;
        fitObj.sse(mwti) = gof.sse;
        
        % formula
        ff = formula(f);
        cn = coeffnames(f);
        cv = coeffvalues(f);
        % replace formala with coefficent 
        for v = 1:numel(cv)
            ff = regexprep(ff,cn{v},num2str(cv(v)));
        end
        forml = ff;
        syms x
        fh(x) = forml
        
        a  = limit(fh,30)
        %%
        

%         for v = 1:numel(cn)
%             syms(cn{v})
%         end
%         symfun(ff,coeffvalues(f))
        
        %%
%         fplot(f)
        return
   end
   
   return
end

return












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
    fprintf('%s %s: %.0f%% significant\n',msr,phenotype,nP*100)
    report(sri,msri) = nP;
end

end

%% export report

a = array2table(report,'VariableNames',msrlist);
cd(pSave); writetable(a,sprintf('Percent significant a95p %s.csv',phenotype));

return






















