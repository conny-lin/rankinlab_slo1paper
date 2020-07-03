function MWTSet = Dance_InitialEtohSensitivityPct(pMWT,varargin)
%% INFORMATION
% chor (added 20170606): built in chor

%% DEFAULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
functionfilepath = mfilename('fullpath');
pSave = '/Users/connylin/Dropbox/RL/Dance Output';
timestamp = 'off';

% settings ++++++
% pvsig = 0.05;
% pvlim = 0.001;
startTime = 90; 
endTime = 95;

Varin = varargin; % store varargin in another variable
vararginProcessor; % run varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% standard Dance processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MWTSet
[MWTSet,pSave] = DanceM_MWTSetStd3(pMWT,functionfilepath,Varin,...
    pSave,'timestamp',timestamp); % standard MWTSET
% chor
pMWTinput = pMWT; % record original pMWT
[~,pMWTpass,pMWTfailed] = chormaster5('TrinityOnly',pMWT,'convert2mat',true); % chor
pMWT = pMWTpass; % update passed pMWT
% update MWTDB 
MWTSet.MWTDB_original = parseMWTinfo(pMWTinput);
MWTSet.MWTDB_failed = parseMWTinfo(pMWTfailed);
MWTSet.MWTDB = parseMWTinfo(pMWT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MWTDB = MWTSet.MWTDB;
DataMeta = IS_getData_Trinity(MWTDB.mwtpath,MWTDB.mwtid, startTime, endTime); % get data from gangnam
IS_stats(DataMeta,MWTDB,pSave); % stats
IS_graph(DataMeta,MWTDB,pSave,startTime,endTime); % calculate initial sensitivity eStats
save(fullfile(pSave,'data.mat'),'DataMeta','MWTDB'); % save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

%% SUB FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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








