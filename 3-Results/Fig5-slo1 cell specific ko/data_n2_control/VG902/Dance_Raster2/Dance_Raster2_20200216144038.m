function MWTSet = Dance_Raster2(pMWT,varargin)

%% setting
% standard dance vairables
timestamp = 'off';
pSave = '/Users/connylin/Dropbox/RL Dance Output';
functionfilepath = mfilename('fullpath');
funName = mfilename;

% savefoldername = mfilename;
graphopt = true; % plot on or off
% create time frames
tap = 1:30;
beforetap = 1;
aftertap = 1;

% standard dance varagin catcher
Varin = varargin;
vararginProcessor; % vararginProcessor 
% ---------------------

%% VALIDATE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DB = parseMWTinfo(pMWT); % get MWT info
% create time frames
% get isi information
ISI = unique(DB.ISI);
if numel(ISI) ~= 1; error('can only accomodate 1 isi'); end
% preplate
preplate = unique(DB.preplate);
if numel(preplate) ~= 1; error('can only accomodate 1 preplate'); end
% tapN
tapN = unique(DB.tapN);
if numel(tapN) ~= 1; error('can only accomodate 1 tapN'); end
if tap(end)>tapN; error('# tap requested is more than actual number of taps');end
% generate the rest
taptimes = (tap*ISI)+preplate-ISI;
startList = taptimes-beforetap; 
endList = taptimes+aftertap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% standard Dance processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MWTSet
[MWTSet,pSave] = DanceM_MWTSetStd3(pMWT,functionfilepath,Varin,...
    pSave,'timestamp',timestamp); % standard MWTSET
% chor
pMWTinput = pMWT; % record original pMWT
[~,pMWTpass,pMWTfailed] = chormaster5('TrinityOnly',pMWT,'deletedat',1); % chor
pMWT = pMWTpass; % update passed pMWT
% update MWTDB 
MWTSet.MWTDB_original = parseMWTinfo(pMWTinput);
MWTSet.MWTDB_failed = parseMWTinfo(pMWTfailed);
MWTSet.MWTDB = parseMWTinfo(pMWT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% process ++++++++++++++++++++++++++
% p = fullfile(pSave,savefoldername);
% pSave = create_savefolder(p);

% get strain +++++++++++++++++++
% a = unique(MWTDB.strain);
% a(regexpcellout(a,'N2')) = [];
% if numel(a) == 1
%     strain = char(a);
% else
%     fprintf('too many strains\n');
%     return 
% end
% -----------------------------------

%% RASTER PLOT - RED + GREEN + ephys data output
MWTDB = MWTSet.MWTDB;
% get unique groups 
gnlist = unique(MWTDB.groupname);
% declare ephys data output
% EphysT = struct;

for gi = 1:numel(gnlist) % cycle through groups

        % prepare local variables
        gn = gnlist{gi}; % get group name
        i = ismember(MWTDB.groupname, gn); % get index to mwt path
        pMWT = MWTDB.mwtpath(i); % get pMWT path
        pSaveG = create_savefolder(pSave,gn); % create group folder

        %% load trinity data ++++++++++++++++
        fprintf('-- loading data...')
        Trinity = cell(numel(pMWT),1);
        Trinity(:,1) = pMWT;
        for mi = 1:numel(pMWT)
            Trn = load_Trinity(pMWT{mi});
%           Trn = get_speeddata_fromTrinity_v3(Trn);
            Trinity{mi,2} = Trn;
        end
        fprintf('finished\n');
        % ------------------------------------

        %% get assay time from trinity data +++++++
% 
%         for ti = 1:numel(startList)
%             start = startList(ti);
%             finish = endList(ti);
% 
%             Data = cell(size(pMWT));
%             mwtid = Data;
%             rTime = mwtid;
% 
%             for mi = 1:numel(pMWT)
%                 %% get plate data
%                 Trn = Trinity{mi,2};            
%                 % get assay time based on tap time
%                 assayTime = trinity_calAssaytime(Trn,start,finish);
%                 % get data
%                 D = trinity_getDataWithinTime(Trn,assayTime);
%                 % translate timeframe
%                 D = trinity_getTimeFrame(D,assayTime,pMWT{mi});
%                 if ~isempty(D)
%                     %% exclude worms
%                     D = exclude_wormNotInAssayTime(D,assayTime);
%                     [D,wormidU] = trinity_avg_speedb(D);
%                     % get time
%                     t = repmat(assayTime(1:end-1),size(D,1),1);
%                     % add to master 
%                     Data{mi} = D;
%                     mwtid{mi} = wormidU;
%                     rTime{mi} = t;
%                 end
%             end
% 
% return
%         end
%         % -----------------------------------------
%         
%         %% get speed data from trinity +++++
%         [SumData,legend] = get_speeddata_fromTrinity_v2(Trinity);
%         
        %% ---------------------------------

        % +++++++++++++++++++++++++++++++++++++++++++++++++++++
        for ti = 1:numel(startList) % cycle through time points 
            % get time
            start = startList(ti); 
            finish = endList(ti); 


            % get data from trinity ++++++++++
            [Data,rTime,savename,mwtid,~] = TrinitySum_multiPlates(pMWT,finish,start,'TrinityData',Trinity);
            % save data
            p = fullfile(pSaveG,[savename,' rasterData.mat']);
            save(p,'pMWT','Data','mwtid','rTime');
            % ----------------------------------

            if graphopt
                % make hot figure
                f1 = rasterPlot_colorSpeed_hot(Data);
                % get print path
                p = fullfile(pSaveG,'raster hot'); if ~isdir(p); mkdir(p); end
                % print
                set(f1,'PaperPositionMode','auto'); % set to save as appeared on screen
                cd(p); print (f1,'-depsc', '-r1200', savename); % save as eps
                print(f1,'-dtiff', savename); % save as (r600 will give better resolution)
                close;


                % make figure (green, jet)
                [f1,Img] = rasterPlot_colorSpeed_gradient2(Data,rTime,0);
                % get print path
                p = fullfile(pSaveG,'raster jet'); if ~isdir(p); mkdir(p); end
                % print
                set(f1,'PaperPositionMode','auto'); % set to save as appeared on screen
                cd(p); print (f1,'-depsc', '-r1200', savename); % save as eps
                print(f1,'-dtiff', savename); % save as (r600 will give better resolution)
                close;


                % make lava plot (sort proportion)
                f1 = moveproportion(Data,rTime,0);
                p = fullfile(pSaveG,'lava proportion plot'); if ~isdir(p); mkdir(p); end
                titlename = sprintf('%s/speed_dir_proportion_%s',p,savename);
                print(f1,'-depsc', titlename); 
                print(f1,'-dtiff', titlename); % save as (r600 will give better resolution)
                close;
            end
%             % electrophys
%             M = nanmean(Data);
%             EphysT(ti).(gn) = M';
%             EphysT(ti).([gn,'_SE']) = (std(Data)./sqrt(size(Data,1)-1))';

        end
%         % write electrophys data
%         p = fullfile(pSaveG,'ephys data.mat');
%         save(p,'EphysT');

end



 



