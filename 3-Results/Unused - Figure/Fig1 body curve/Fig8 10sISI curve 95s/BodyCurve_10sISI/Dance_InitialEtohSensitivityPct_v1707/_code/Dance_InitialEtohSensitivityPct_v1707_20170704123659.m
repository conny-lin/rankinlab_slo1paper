function MWTSet = Dance_InitialEtohSensitivityPct_v1707(pMWT,pSave,varargin)
%% INFORMATION
% chor (added 20170606): built in chor

%% DEFAULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functionfilepath = mfilename('fullpath'); % function paths
% funname = mfilename; % funciton name
% timestamp = 'off';
% pSave = '/Users/connylin/Dropbox/RL/Dance Output';
% settings ++++++
% pvsig = 0.05;
% pvlim = 0.001;
startTime = 90; 
endTime = 95;

V = varargin; % store varargin in another variable
vararginProcessor; % run varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% standard Dance processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MWTSet,pSave] = DanceM_MWTSetStd_v1707(pMWT,mfilename('fullpath'),V,pSave); % standard MWTSET
[Legend,pMWTpass,pMWTfailed] = chormaster5('TrinityOnly',pMWT,'convert2mat',true); % chor
[MWTSet,pMWT] = DanceM_MWTSet_postChor(MWTSet,pMWT,pMWTpass,pMWTfailed,Legend); % update MWTSet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting -----------------------------------------------------------------
MWTDB = MWTSet.MWTDB;
%--------------------------------------------------------------------------
% get data ----------------------------------------------------------------
DataMeta = IS_getData_Trinity_v1707(MWTDB, startTime, endTime); % get data from gangnam
IS_stats(DataMeta,MWTDB,pSave); % stats
IS_graph(DataMeta,MWTDB,pSave,startTime,endTime); % calculate initial sensitivity eStats
%--------------------------------------------------------------------------

%% PCT STATS 

CS = CurveStats_v1707;
CS.Data = DataMeta;
CS.MWTDB = MWTDB;

fid = fopen(fullfile(pSave,'curve pct.txt'),'w');
fprintf(fid,'%s',CS.statstr);
fclose(fid);

MWTSet.Stats = CS;
% SAVE 
save(fullfile(pSave,'data.mat'),'MWTSet','DataMeta','MWTDB'); % save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

%% SUB FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function DataMeta = IS_getData_Trinity(pMWT,mwtID,startTime,endTime)
% 
%     load('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/Chor_output_handling/legend_trinity.mat')
%     legend_trinity(end+1) = {'id'};
%     DataMeta = cell(size(pMWT)); % create array
%     mwtIDheader = cell(size(pMWT)); % create mwtid header
% 
% 
%     for mwti = 1:numel(pMWT)
%         processIntervalReporter(numel(pMWT),20,'MWT',mwti);
%         pmat = char(getpath2chorfile(pMWT(mwti),'trinitySummary.mat','reporting',0)); % get pmat
%         D = load(pmat); % load data
% 
%         % extract time
%         D1 = D.masterData(:,2);
%         nWorms = size(D1,1);
%         ti = nan(nWorms,1);
%         tf = ti;
%         for wrmid = 1:nWorms
%             ti(wrmid) = D1{wrmid,1}(1,1);
%             tf(wrmid) = D1{wrmid,1}(end,1);
%         end
%         D = D.masterData(ti <= startTime & tf >= endTime,:); % get data within time frame
% 
%         if ~isempty(D) 
% 
%             % add worm id to last column
%             ci = strcmp(legend_trinity,'id');
%             a = D(:,1);
%             b = regexprep(a,'^0{1,}','');
%             id = cellfun(@str2num,b); % worm id
%             nWorms = size(D,1);        
%             for wrmi = 1:nWorms
%                 id1 = id(wrmi);
%                 n = size(D{wrmi,2},1);
%                 D{wrmi,2}(:,ci) = repmat(id1,n,1);
%             end 
% 
%             % convert other datadata
%             D = cell2mat(D(:,2)); 
%             t = D(:,strcmp(legend_trinity,'time'));
%             D(~(t >= startTime & t <= endTime),:) = [];
% 
%             % calculate speed, speedbm, curve by worms
%             speed = D(:,strcmp(legend_trinity,'speed'));
%             id = D(:,strcmp(legend_trinity,'id'));
%             midline = D(:,strcmp(legend_trinity,'midline'));
%             curve = D(:,strcmp(legend_trinity,'curve'));
%             speedbm = speed./midline;
% 
%             clear S
%             S = grpstatsTable(speed,id,'gnameutitle','wrmid');
%             S.wrmid = cellfun(@str2num,S.wrmid);
%             SB = grpstatsTable(speedbm,id,'gnameutitle','wrmid');
%             SB.wrmid = cellfun(@str2num,SB.wrmid);
%             C = grpstatsTable(curve,id,'gnameutitle','wrmid');
%             C.wrmid = cellfun(@str2num,C.wrmid);
% 
%             if ~isequal(S.wrmid, SB.wrmid) || ~isequal(SB.wrmid,C.wrmid)
%                 error('wormid does not match');
%             end
%             D1 = [S.wrmid C.mean SB.mean S.mean];
% 
%             % put in main datasheet
%             DataMeta{mwti} = D1;
%             mwtIDheader{mwti} = repmat(mwtID(mwti),size(D1,1),1);
% 
%         end
%     end
%     % collapse data
%     Data = [cell2mat(mwtIDheader) cell2mat(DataMeta)];
%     DataMeta = array2table(Data,'VariableNames',{'mwtid','wrmid','curve','speedbm','speed'});
% end


% function IS_stats(DataMeta,MWTDB,pSave)
%     gns = MWTDB.groupname(DataMeta.mwtid);
% 
%     msrlist = {'curve','speedbm','speed'};
%     for msri = 1:numel(msrlist)
%         msr = msrlist{msri};
%         [strain,rx] = parse_groupname(gns);
%         rx(cellfun(@isempty,rx)) = {'0mM'};
%         G = [{strain} {rx}];
%         gvar = {'strain','dose'};
%         Y = DataMeta.(msr);
%         anovan_std(Y,G,gvar,pSave,msr)
%     end
% end


% function IS_graph(DataMeta,MWTDB,pSave,startTime,endTime)
% gns = MWTDB.groupname(DataMeta.mwtid);
% gnseq = unique(gns);
% i = regexpcellout(gnseq,'N2');
% gnseq = [gnseq(i);gnseq(~i)];
% MWTSet.Data_GroupByWorm = struct;
% MWTSet.Data_GroupByWorm.curve = grpstatsTable(DataMeta.curve, gns,'gnameu',gnseq,'gnameutitle','groupname');
% MWTSet.Data_GroupByWorm.speedbm = grpstatsTable(DataMeta.speedbm, gns,'gnameu',gnseq,'gnameutitle','groupname');
% MWTSet.Data_GroupByWorm.Speed = grpstatsTable(DataMeta.speed, gns,'gnameu',gnseq,'gnameutitle','groupname');
% 
% %% graph
% GraphData = MWTSet.Data_GroupByWorm;
% msrlist = fieldnames(MWTSet.Data_GroupByWorm);
% for msri = 1:numel(msrlist)
%     close;
%     gn = GraphData.(msrlist{msri}).groupname;
%     gn = regexprep(regexprep(gn,'_',' '),'mM','');
%     y = GraphData.(msrlist{msri}).mean;
%     e = GraphData.(msrlist{msri}).se;
%     msr = msrlist{msri};
%     savename = sprintf('%s t%d-%d',msr,startTime,endTime);
% 
%     % plot
%     fig1 = figure('Visible','off');
%     ax1 = axes('Parent',fig1,'Box','off','XTick',1:numel(gn),...
%         'XTickLabel',gn');
%     hold(ax1,'on');
% 
%     bar(y,'EdgeColor','none','FaceColor',[0 0.447058826684952 0.74117648601532])
%     hold on;
%     errorbar([1:numel(gn)]',y,e,'Color',[0 0 0],'LineStyle','none','LineWidth',1.2)
% 
%     ylabel(msrlist{msri})
%     xlim([.5 numel(gn)+.5])
%     printfig(savename,pSave,'w',3.5,'h',2,'closefig',1)
% end
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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








