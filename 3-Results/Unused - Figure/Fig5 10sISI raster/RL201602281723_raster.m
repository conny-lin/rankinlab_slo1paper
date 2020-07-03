%% get raster plot for N2 valid exp
%% INITIALIZING
clc; clear; close all;
%% SAVE PATHS
pSave = mfilename('fullpath'); if isdir(pSave) == 0; mkdir(pSave); end
%% ADD FUNCTION PATH
pFun = {'/Users/connylin/Dropbox/Code/Matlab/Library/General';
        '/Users/connylin/Dropbox/Code/Matlab/Library';
        '/Users/connylin/Dropbox/RL/Code/Modules'};
for x = 1:numel(pFun); 
    p = pFun{x};
    if isdir(p) == 0; error('function path %d not found'); end
    addpath(p); 
end
%% SETTING
% pData = '/Volumes/COBOLT/MWT';
% pD = /Dance_Glee_Showmance';
pSaveHome = '/Users/connylin/Dropbox/RL/PhD/Chapters/2-STH N2/Data/10sISI/2-Valid Exp';
NWORMS = Inf;
expLimit = 'N2 within';
frameInt = 10;

% create time frames
startList = [95:10:400];
assaywindow = 10;
endList = startList+assaywindow;

%% LOAD MWT DATABASE
load('/Users/connylin/Dropbox/rl/Publication/PhD Dissertation/Chapters/3-STH N2/3-Results/10sISI 400mM/m1-sample/select_valid_sample_10sISI_RL201602201421/MWTInfo.mat');
pMWT = MWTDB.mwtpath;

%% CHOR
pMWTc = convertTrinityDat2Mat(pMWT,1); 
L = chormaster4('Trinity',pMWTc);
% summarize trinity data and delete .dat file to save memory
pMWTbad = convertTrinityDat2Mat(pMWTc,1); 
% exclude bad files
pMWToriginal = pMWT;
pMWT(ismember(pMWT,pMWTbad)) = [];
Db = parseMWTinfo(pMWT);
groupnameList = unique(Db.groupname);

%% raster plot per time points
pSaveA = pSave;
for ti = 18:numel(startList)
    start = startList(ti);
    finish = endList(ti);        
    T = table; % for electrophys
    for gi = 1:numel(groupnameList)
        gT = groupnameList{gi};
        fprintf('processing: %s %d-%ds\n',gT,start,finish)
        pMWTG = pMWT(ismember(Db.groupname,gT));
        pSaveG = [pSaveA,'/',gT];
        if isdir(pSaveG) == 0; mkdir(pSaveG); end
        % run raster
        [f1,Data,savename,Tog,mwtid,rTime,~] = rasterPlot_colorSpeed(pMWTG,start,finish,...
            'NWORMS',NWORMS,'visibleG',0);
        % save fig
        cd(pSaveG);
        set(f1,'PaperPositionMode','auto'); % set to save as appeared on screen
        print (f1,'-depsc', '-r1200', savename); % save as eps
        close;
        % save data
        save(sprintf('%s/%s rasterData.mat',pSaveG,savename),'pMWT','Data','mwtid','Tog','rTime');
   
        % electrophys
        D = Data;
        M = mean(D);
        T.(gT) = M';
        T.([gT,'_SE']) = (std(D)./sqrt(size(D,1)-1))';
    end
    writetable(T,sprintf('%s/ef graph %d-%d.csv',pSaveA,start,finish));
end
%% done
fprintf('\n\nDONE\n'); return














