%% SAVE PATHS
pSave = mfilename('fullpath'); if isdir(pSave) == 0; mkdir(pSave); end
%% ADD FUNCTION PATH
pFun = {'/Users/connylin/Dropbox/Code/Matlab/Library/General';
        '/Users/connylin/Dropbox/Code/Matlab/Library/Graphs';
        '/Users/connylin/Dropbox/rl/Code/Modules/MWTDatabase';
        '/Users/connylin/Dropbox/rl/Code/Modules/Graphs/starfish';
        '/Users/connylin/Dropbox/RL/Code/Modules/Chor'};
for x = 1:numel(pFun); 
    p = pFun{x};
    if isdir(p) == 0; error('function path %d not found'); end
    addpath(p); 
end


%% input variables
timeAssay = [99:10:98+(30*10)];
timeAssay(2,:) = timeAssay+ 3;

assayInt = 0.2;
Legend = chormaster4('Gangnam',{});
Legend =  regexprep(Legend{1}{2},':','_');
nwormselect = 100;

%% Load sample pMWT paths
pdata = '/Users/connylin/Dropbox/RL/Publication/PhD Dissertation/Chapters/3-STH N2/2-Method/worm path/pMWT.mat';
load(pdata,'pMWT');


% group pMWT
MWTDB = parseMWTinfo(pMWT);
gnameu = unique(MWTDB.groupname);
pMat = cellfun(@strcat,pMWT,cellfunexpr(pMWT,'/Gangnam.mat'),'UniformOutput',0);
pMatG = cell(size(gnameu));
for gi = 1:numel(gnameu)
    pMatG{gi} = pMat(ismember(MWTDB.groupname,gnameu{gi}));
end


%% group trinity data
        
% get trinity files within specific time frame
GroupData = cell(size(gnameu));
for gi = 1:numel(gnameu)
    gn = gnameu{gi};
    for ti = 1:size(timeAssay,2)
        % time variables
        assayTime1 = timeAssay(1,ti);
        assayTime2 = timeAssay(2,ti);
        assayTimeseries = assayTime1:assayInt:assayTime2;
        fprintf('processing: %s %.0fs-%.0fs\n',gn,assayTime1,assayTime2);
        
        pMat = pMatG{gi};
        nFiles = numel(pMat);
        XG = cell(nFiles,1);
        BiasG = cell(nFiles,1);
        TapG = cell(nFiles,1);
        YG = cell(nFiles,1);
        wormidG = cell(nFiles,1);
        for fi =1:numel(pMat)
            % process each trinity files
            pfile = pMat{fi};
            Import = load(pfile);

            % get only time wtihin assay time
            iworm = Import.time(:,1) <=assayTime1 & Import.time(:,2) >=assayTime2;
            Data = array2table(cell2mat(Import.Data(iworm)),'VariableNames',Legend);
            % get rid of unwanted times
            Data(Data.time < assayTime1 | Data.time > assayTime2,: ) = [];
            
            % input data
            plateid = fi*100000000;
            wormidG{fi} = Data.id+plateid;
            XG{fi} = Data.loc_x;
            YG{fi} = Data.loc_y;
            BiasG{fi} = Data.bias;
            TapG{fi} = Data.tap;
        end
        % convert to numeric array
        wormid = cell2mat(wormidG);
        X = cell2mat(XG);
        Tap = cell2mat(TapG);
        Y = cell2mat(YG);     
        Bias = cell2mat(BiasG);
        
        % randomly chose 100 worms
        wormidu =  unique(wormid);
        nworm = numel(wormidu);
        if nworm < nwormselect
            error('worm N smaller than wanted');
        end
        nselect = randi(nworm,nwormselect,1);
        wormidselect = wormidu(nselect);
        iselect = ismember(wormid,wormidu(nselect));
        wormid = wormid(iselect);
        X = X(iselect);
        Y = Y(iselect);
        Bias = Bias(iselect);
        Tap = Tap(iselect);
        
        % figure
        fig1 = starFish(wormid,X,Y,Bias,Tap,'visibleopt',0,'tapcirclesize',10,'pathlinewidth',1,'wormidshow',0,'alignAtTap',1);
        % label and save
        tname = sprintf('%s %.1fs-%.1fs N(worm)=%d',regexprep(gn,'_','-'),assayTime1,assayTime2,nwormselect);
        title(tname)
        savename = sprintf('wormpath %s %.0fs-%.0fs random N_%d',regexprep(gn,'_','-'),assayTime1,assayTime2,nwormselect);
        savefigepsOnly(savename,pSave);
    end
end




fprintf('\nDone\n');


























