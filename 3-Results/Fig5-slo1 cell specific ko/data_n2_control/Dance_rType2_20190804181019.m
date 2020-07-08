function MWTSet = Dance_rType2(pMWT,varargin)

%% varargin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard dance vairables
timestamp = 'off';
pSave = '/Users/connylin/Dropbox/RL Dance Output';
functionfilepath = mfilename('fullpath');
% time_tap_col = 6;
time_baseline_col = 3:5;
time_response_col = 7:10;
rTarget = 'accelerate forward';
rTargetName = 'acceleration';
genotype = '';
strain = '';
n_lowest = 1;
msrlist = {'AccProb'};
%% standard dance varagin catcher
Varin = varargin;
vararginProcessor; % vararginProcessor 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% standard Dance processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MWTSet
[MWTSet,pSave] = DanceM_MWTSetStd3(pMWT,functionfilepath,Varin,...
    pSave,'timestamp',timestamp); % standard MWTSET
% chor
% pMWTinput = pMWT; % record original pMWT
% [~,pMWTpass,pMWTfailed] = chormaster5('TrinityOnly',pMWT,'deletedat',1); % chor
% pMWT = pMWTpass; % update passed pMWT
% % update MWTDB 
% MWTSet.MWTDB_original = parseMWTinfo(pMWTinput);
% MWTSet.MWTDB_failed = parseMWTinfo(pMWTfailed);
% MWTSet.MWTDB = parseMWTinfo(pMWT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  Dance Raster     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pSave_Raster = fileparts(pSave); % create save folder
Data_Raster = Dance_Raster2(pMWT,'pSave',pSave_Raster,'graphopt',false); % Dance raster
MWTSet.PATHS.pSave_Raster = Data_Raster.PATHS.pSaveA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pData = MWTSet.PATHS.pSave_Raster;
[~,~,groupnames] = dircontent(pData);

MeanMaster = cell(numel(groupnames),1);
NMaster = MeanMaster;
TimeMaster = MeanMaster;
pMWTMaster = MeanMaster;

for gi =1:numel(groupnames)

    gn = groupnames{gi};
    p = fullfile(pData,gn);
    rasterfiles = dircontent(p,'rasterPlot*mat');

    % find t1
    a = regexpcellout(rasterfiles,'(_)|\s','split');
    t1 = cell2mat(cellfun(@str2num,a(:,2),'UniformOutput',0));
    tf = cell2mat(cellfun(@str2num,a(:,3),'UniformOutput',0));
    n = cell2mat(cellfun(@str2num,a(:,5),'UniformOutput',0));

    [c,seq] = sort(t1);
    rasterfiles = rasterfiles(seq);

    MeanSum = cell(numel(rasterfiles),1);
    SESum = MeanSum;
    NSum = MeanSum;
    MWTPATH = MeanSum;
    TimeSum = MeanSum;

    for ri = 1:numel(rasterfiles)
        % get data -------------------------
        p = fullfile(pData, gn,rasterfiles{ri});
        Data = load(p);
        D = Data.Data;
        rTime = Data.rTime(1,:);
        mwtid = Data.mwtid;
        pMWT = Data.pMWT;
        MWTDB = parseMWTinfo(pMWT);
        % -----------------------------------

        % get response type summary %--------------
        Baseline = D(:,time_baseline_col);
        Response = D(:,time_response_col);
        [T1, RTclean, R, RN,leg] = compute_response_type(Response, Baseline);
        R = cell2table(R);
        % ---------------------------------

        % add plate info % --------------
        plateinfo = table;
        plateinfo.mwtid = mwtid;
        plateinfo.timeid = repmat(ri,size(mwtid,1),1);
        plateinfo.wormid = (1:numel(mwtid))';
        % ---------------------------------

        % calculate response type (acceleration) --------------
        % summarize response type by plate
        [A,N,~,mwtidlist] = responseType_pct_plate(plateinfo,R,rTarget,'Type','any');
        % rid of data with NaN ouputs
        invalid = any(N < n_lowest,2) & any(isnan(A),2);
        % rid of N lower than n limist
        RespNv = A;
        RespNv(invalid,:) = NaN;
        % ---------------------------------------

        MWTPATH{ri} = pMWT(mwtidlist);
        MeanSum{ri} = RespNv;
        NSum{ri} = N;
        TimeSum{ri} = repmat(ri,size(A,1),1);
    end

    % take out cell
    MeanSum = cell2mat(MeanSum);
    NSum = cell2mat(NSum);
    TimeSum = cell2mat(TimeSum);
    MWTPATH = celltakeout(MWTPATH);

    MeanMaster{gi} = MeanSum;
    NMaster{gi} = NSum;
    TimeMaster{gi} = TimeSum;
    pMWTMaster{gi} = MWTPATH;



end


MeanMaster = cell2mat(MeanMaster);
NMaster = cell2mat(NMaster);
TimeMaster = cell2mat(TimeMaster);
pMWTMaster = celltakeout(pMWTMaster);
% =========================================================================


%% PREPARE DATA TABLE =====================================================
M = parseMWTinfo(pMWTMaster);
D = table;
D.mwtpath = M.mwtpath;
D.groupname = M.groupname;
D.strain = M.strain;
% create variables
rx = M.rx;
ictrl = strcmp(rx,'NA');
rxe = rx(~ictrl);
if sum(regexpcellout(rxe,'mM')) == numel(rxe)
    D.dose = rx;
    D.dose(ismember(D.dose,'NA')) = {'0mM'};
else
    D.rx = rx;
end
% examine variable types
factors = {};
if ismember('dose',D.Properties.VariableNames)
    factors{end+1} = 'dose';
else
    factors{end+1} = 'rx';
end
if numel(unique(D.strain)) >1
    factors{end+1} = 'strain';
end

% dose = a(:,2);
% dose(cellfun(@isempty,dose)) = {'0mM'};
% D.dose = dose;
% D.strain = strain;

D.tap = TimeMaster;
D.AccProb = MeanMaster;
D.n = NMaster;
Data = D;

clear D;
% =========================================================================


%% graph ===============================
G = statsByGroup(Data, msrlist,'tap','groupname'); % get graphic data 
for msri = 1:numel(msrlist)

    msr = msrlist{msri};
    
    % get graph data ----------------
    D = G.(msr);
    y = D.Y;
    e = D.E;
    x = D.X;
    gn = D.groupname;
    D1(:,:,1) = x;
    D1(:,:,2) = y;
    D1(:,:,3) = e;
    % sort by N2
    gns = sortN2first(gn,gn);
    [i,j] = ismember(gn,gns);
    D1 = D1(:,j,:);
    X = D1(:,:,1);
    Y = D1(:,:,2);
    E = D1(:,:,3);
    % -------------------------------

    % setting --------------------------
    w = 4.5;
    h = 3.5;
    gp = graphsetpack('cathyline');
    gnss = regexprep(gns','_',' ');
    gnss = regexprep(gnss,strain,genotype);
    gnss = regexprep(gnss,'N2','Wildtype');
    gp.DisplayName = gnss;
    gpn = fieldnames(gp);
    % -----------------------------------

    % plot ------------------------------------------
    fig = figure('Visible','off','PaperSize',[w h],'Unit','Inches',...
        'Position',[0 0 w h]); 
    ax1 = axes('Parent',fig,'XTick',[0:10:30]); 
    hold(ax1,'on');
    e1 = errorbar(X,Y,E,'Marker','o','MarkerSize',4.5);
    % -----------------------------------------------

    % get settings -----------------------------------
    for gi = 1:numel(gn)
        for gpi = 1:numel(gpn)
            e1(gi).(gpn{gpi}) = gp.(gpn{gpi}){gi};
        end
    end
    % -----------------------------------------

    % title -------------------------
    title(genotype)
    % -------------------------------

    % legend ------------------------
    %     lg = legend('show');
    %     lg.Location = 'eastoutside';
    %     lg.Box = 'off';
    %     lg.Color = 'none';
    %     lg.EdgeColor = 'none';
    %     lg.FontSize = 8;
    % -------------------------------

    % x axis -----------------------------------------------
    xlabel('Tap')
    xlim([0 30.5]);
    % --------------------------
    % y axis ---------------------
    yname = sprintf('P (%s)',rTargetName);
    ylim([0 1]);
    ylabel(yname);
    % ------------------------
    % axis ------------------------------------
    axs = get(ax1);
    ax1.TickLength = [0 0];
    ax1.FontSize = 10;
    % ----------------------------------------------------

    % save -------------------------------------------
    savename = sprintf('%s/%s',pSave,msr);
    printfig(savename,pSave,'w',w,'h',h,'closefig',1);
    % write in csv
    T = saveGraphstruct2table(G.(msr),fullefile(pSave,[msr,'.csv'])); %added 20190804
    % -------------------------------------------
end



% =========================================================================



%% DESCRIPTIVE STATS & RMANOVA N=PLATES ==================================
% anova settings +++++++++++++++++
MWTDB = MWTSet.MWTDB;
rmName = 'tap';
idname = 'mwtpath';
% fNameInd = 'groupname';
% fNameS = strjoin(factors,'*');


% descriptive stats +++++++++++++++++++++++++++++++
txt1 = '*** Descriptive Stats ***';
t = tabulate(MWTDB.groupname);
% gname
gname = sortN2first(t(:,1),t(:,1));
gname = strjoin(gname,', ');
txt2 = sprintf('gname = %s',gname);
% n
n = sortN2first(t(:,1),t(:,2));
n = strjoin(cellfun(@num2cellstr,n),', ');
txt3 = sprintf('N(plates) = %s',n);
Ntext = sprintf('%s\n%s\n%s',txt1,txt2,txt3);
% -----------------------------------------------

for msri = 1:numel(msrlist)

    msr = msrlist{msri};

    % anova +++++++++++++++++++++++++++++++++++++++
%     A = anovarm_convtData2Input(Data,rmName,[fName,fNameInd],msr,idname);
%     textanova = anovarm_std(A,rmName,fNameInd,fNameS);
  [textanova,DS] = anovarm_std(Data,rmName,factors,idname,msr);
    % -----------------------------------------------

    textout = sprintf('%s\n%s',Ntext,textanova); % join
    
    % save data ++++++++++++
    cd(pSave);
    p = sprintf('%s RMANOVA.txt',msr);
    fid = fopen(p,'w');
    fprintf(fid,'%s',textout);
    fclose(fid);
    % ----------------------

end
% =========================================================================


%% save data ++++++++++++
cd(pSave); save('data.mat','Data','G','MWTDB','MWTSet');
% ----------------------

















