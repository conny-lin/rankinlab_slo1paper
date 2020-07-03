
%% INITIALIZING
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
addpath('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/Dance/Dance_RespType');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% ---------------------------

%% GLOBAL INFORMATION
% paths & settings -----------
pSave = fileparts(pM);
addpath(fileparts(pM));
pData = '/Users/connylin/Dropbox/RL/RL Pub PhD Dissertation/Chapters/2-Wildtype/3-Results/1-Figures & Data/Fig2-10 raster 10sISI/raster per tap ephys_201602281723';
% '/Users/connylin/Dropbox/RL Pub PhD Dissertation/Chapters/4-EARS/3-Results/1-N2 EARS/Data/raster per tap ephys_201602281723';
% ---------------------------

% strains %------
% get group list
[~,~, glist] = dircontent(pData);
glist(regexpcellout(glist,'graffle')) = [];
% strainlist = dircontent(pData);
% get strain info
% strainNames = DanceM_load_strainInfo(strainlist);
%----------------

% settings %--------------------
n_lowest = 10;
% time_baseline = [-0.3 -0.1];
tapstarttime = 100;
isi = 10;
tapN = 30;
assaytimeperiod = 1; % 1s after tap
baselinetimeperiod = 0.3; % 0.3 seconds before tap
responseTypes = {'accelerate forward'; 'reversal';'pause'};
rTarget = {'accelerate forward'; 'reversal'};
% construct times
taptimelist = [tapstarttime:isi:(tapstarttime+isi*(tapN-1))];
matfiletimelist = taptimelist-5;
% time_response = [0.1 0.5];
% statsName = 'exclude no response';
% 
% % create time frames
% startList = [98 388];
% endList = startList+assaywindow;
% --------------------------------



%% RESPONSE PROBABILITY: ACC, PAUSE, NO RESPONSE, REVERSAL, DECELLEARTION
% construct output array -------------
timecolname = matfiletimelist';
ncol = numel(matfiletimelist);
nrow = numel(glist);
Out = struct;
n = numel(matfiletimelist)*numel(responseTypes)*numel(numel(glist));
DataStats = cell(n,1);
iDataStats =1;
% -----------------------------

for gi = 1:numel(glist) % cycle through groups
    
    % report progress %------------
    fprintf(''); % separator
    processIntervalReporter(numel(glist),1,'group',gi);
    % ------------------------------
        
    % get .mat file list and time list --------------
    groupname = glist{gi};
    matfilelist = dircontent(fullfile(pData, groupname),'*.mat');
    % create mat file table
    matfileinfo = table;
    matfileinfo.fname = matfilelist;
    a = regexpcellout(matfilelist, '_','split');
    matfileinfo.t1 = cellfun(@str2num,a(:,2));
    matfileinfo.tf = cellfun(@str2num,a(:,3));
    n = a(:,5);
    matfileinfo.n = cellfun(@str2num,regexpcellout(n,'^\d{1,}','match'));
    % ----------------------------------------------
    
    % declare array ---------------------
    MEAN = nan(numel(responseTypes), ncol);
    SEM = MEAN;
    NVAL = MEAN;
    % -----------------------------------
    
    for ti = 1:numel(matfiletimelist) % cycle through time choices

        % report progress %------------
        fprintf(''); % separator
        processIntervalReporter(numel(matfiletimelist),5,'-- time',ti);
        % ------------------------------

        % time info -----------------
        matfiletime1now = matfiletimelist(ti);
        timestartnow = taptimelist(ti);
        timeendnow = timestartnow + assaytimeperiod;
        % --------------------------

        % data % -------------
        f = matfileinfo.fname{matfileinfo.t1 == matfiletime1now};
        p = fullfile(pData, groupname,f);
        load(p,'rTime','Data','pMWT','mwtid'); % load data
        rTime = rTime(1,:); % simply rTime
        % find tap time column
        assaytimei = find(rTime >= timestartnow & rTime <= timeendnow);
        responsei = assaytimei(1);
        % get data before tap Baseline
        baselinetimei = responsei-3:1:responsei-1;
        baselinetime = rTime(baselinetimei);
        Baseline = Data(:,baselinetimei);
        % get response data
        responsetimei = responsei:1:responsei+5;
        responsetime = rTime(responsetimei);
        Response = Data(:,responsetimei); 
        % ----------------

        % get response type summary %--------------
        [~,~,R,~,leg] = compute_response_type(Response, Baseline);
        R = cell2table(R);
        % ---------------------------------

        % add plate info % --------------
        MWTDB = parseMWTinfo(pMWT);
        plateinfo = table;
        plateinfo.mwtid = mwtid;
        plateinfo.mwtpath = pMWT(mwtid);
        % ---------------------------------
        
        
        for ri = 1:numel(responseTypes)
            rn = responseTypes{ri};
            % calculate response type (acceleration) --------------
            % summarize response type by plate
            [A,N,NV,mwtidlist] = responseType_pct_plate(plateinfo,R,rn);
            % rid of data with NaN ouputs
            invalid = any(N < n_lowest,2) & any(isnan(A),2);
            % rid of N lower than n limist
            RespNv = A;
            RespNv(invalid,:) = NaN;      
            
            %% get only times of interest (first 0.2s)
            t = timestartnow:0.2:(timestartnow+(5*0.2)); % get times
            t = t(1);
            RespNv = RespNv(:,1); % get only first column of response
            
            % transform data
            a = reshape(RespNv,numel(RespNv),1);
            times = repmat(t,size(RespNv,1),1);
            rt = repmat(ri,size(RespNv,1),1);            
            g = repmat(gi,size(RespNv,1),1);
            DataStats{iDataStats} = [rt g mwtidlist times a];
            iDataStats = iDataStats+1;
            
            % descriptive stats 
            [d,leg] = statsBasic(RespNv);
            m = d(:,3);
            se = d(:,5);
            n = d(:,1);
            
            % place data in output array
            [i,j] = ismember(t,taptimelist);
            MEAN(ri,j) = m;
            SEM(ri,j) = se;
            NVAL(ri,j) = n;
            % ------------------------------------
            
        end 
    end
    

    
    %% organize descripive output data
    Out.(groupname).n = NVAL;
    Out.(groupname).mean = MEAN;
    Out.(groupname).se = SEM;
    
    % table output for descriptive data
    T1 = statsTable_multvar(NVAL,MEAN,SEM,responseTypes);
    T = table;
    T.tap = (1:numel(taptimelist))';
    T = [T T1];
    writetable(T,fullfile(pM,sprintf('%s.csv',glist{gi})));
end


% organize anova stats
DataStats = reshape(DataStats,numel(DataStats),1);
d = cell2mat(DataStats);
T = array2table(d); 
T.Properties.VariableNames = {'response_type', 'group','mwtid','tap','probability'};
T.response_type = responseTypes(T.response_type);
[i,j] = ismember(T.tap,taptimelist);
T.tap = j;
T.group = glist(T.group);  
DataStats = T;

% save
save(fullfile(pM,'data.mat'),'Out','DataStats','responseTypes');



%% FIGURES
% settings %--------------------
% timelist = [95:10:29*10+95];
tapnum = 1:numel(taptimelist);
% --------------------------------

% extract only 0.1s data
glist = fieldnames(Out);
for gi = 1:numel(glist)
    gn = glist{gi};
end
% --------------------------------

% plot errorbar across time (manual)
titlename = regexprep(regexprep(glist,'N2','wild-type'),'_',' ');
[i,j] = ismember(responseTypes,rTarget);
iTarget = j(i);

for gi = 1:numel(glist)
    gn = glist{gi}; % get group specific info 
    d = Out.(gn); % get data 
    % get variables -----------------------------
    y = d.mean(iTarget,:)';
    e = d.se(iTarget,:)';
    x = repmat(tapnum,size(d.mean,1)-1,1)';
    % ------------------------------------------
    fig_pResp(pM,x,y,e,titlename{gi},gn,'dispname',rTarget); % save figure
end



%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% textout = anovarm_std(DataStats,'tap',{'response_type','group'},'probability')
% function textout = anovarm_std(DataTable,rmName,factors,varargin)
DataTable= DataStats;
DataTable(ismember(DataTable.response_type,rTarget),:) = []; % delete irrelevant data

rmName = 'tap';
factors = {'group','response_type'};
dvname = 'probability';
factorcombo = strjoin(factors,' * ');
% varargin ++++++++++++++++++++++++++++++
alpha = 0.05;
pvlimit = 0.001;

% create rmanova options =================================================
rmlist = unique(DataTable.(rmName));
rmvartable = array2table(rmlist,'VariableNames',{rmName});
% ----------------------------------------------


% deal with multiple factors
if numel(factors)>1
    factorName = char(strjoin(factors,'*'));
    g1 = strjoinrows(DataTable(:,factors),'*');
else
    g1 = DataTable.(factors);
end
DataTable.factors = g1;
g = unique(g1);
gg = pairwisecomp_getpairs(g);
gpairs = strjoinrows(gg,' x ');
% -------------------------------------------

%% transpose %=============================================================
% create rows
var = [factors {idname}];
IV = DataTable(:,var);
for i = 1:numel(var)
   if isnumeric(IV.(var{i}))
      IV.(var{i}) = num2cellstr(IV.(var{i}));
   end
end
iv = strjoinrows(IV,'*');
ivu = unique(iv); % get uniqiue dv pairs
% break up ivu
ivusep = regexpcellout(ivu,'*','split');
IVT = cell2table(ivusep,'VariableNames',var);
IVT.factors = ivu;
% assign rows
[~,rowi] = ismember(iv,ivu);

% create output matrix
O = nan(numel(ivu),numel(rmlist));
ind = sub2ind(size(O),rowi,DataTable.(rmName));
O(ind) = DataTable.(dvname);

% create variable names
a = repmat({rmName},numel(rmlist),1);
b = num2cellstr(rmlist);
v = strjoinrows([a b],'');
OT = array2table(O,'VariableNames',v);

DatarmANOVA = [IVT OT];
%=========================================================================


%% RMANOVA ================================================================
D = DatarmANOVA;
% decide which to run
% if multiStrain && multiDose; factorName = 'strain*dose';
% elseif multiStrain && ~multiDose; factorName = 'strain';
% elseif ~multiStrain && multiDose; factorName = 'dose';
% end


% rmanova multi-factor +++++++
rmTerms = sprintf('%s%d-%s%d ~ %s',rmName,rmlist(1),rmName,rmlist(end),factorcombo);
rm = fitrm(D,rmTerms,'WithinDesign',rmvartable);
t = anovan_textresult(ranova(rm), 0, 'pvlimit',pvlimit);
rmText = sprintf('*** RMANOVA(%s:%s) ***\n%s\n',rmName,factorName,t);
% ----------------------------

%% rmanova pairwise by each factor ++++++++
textout ='';
for i = 1:numel(factors)
    compName = factors{i};
%     rmTerms = sprintf('%s%d-%s%d~%s',rmName,rmlist(1),rmName,rmlist(end),compName);
%     rm = fitrm(D,rmTerms,'WithinDesign',rmvartable);
    t = multcompare(rm,compName);
    % text output
    textout = sprintf('%s\n*** Posthoc(Tukey) by %s ***\n',textout,compName);
    if isempty(t)
        textout = sprintf('%sAll comparison = n.s.\n',textout);
    else
        textout = sprintf('%s%s\n',textout,...
            multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha));
    end
end
rmByFactortxt = textout;
% ----------------------------------

%% comparison by taps by group +++++++++++
t = multcompare(rm,compName,'By',rmName);
%  keep only unique comparisons
a = strjoinrows([t.([compName,'_1']) t.([compName,'_2'])],' x ');
t(~ismember(a,gpairs),:) =[]; 
% record
    textout = sprintf('%s\n\n*** Posthoc(Tukey)%s by %s ***\n',textout,rmName,compName); 
if isempty(t)
    textout = sprintf('%s\nAll comparison = n.s.\n',textout);
else
    t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alpha);
    textout = sprintf('%s\n%s\n',textout,t);
end
% --------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% report done --------------
fprintf('\n--- Done ---\n\n' );
return

%% --------------------------























%% plot errorbar across time (manual)
titlename = {'Wildtype 0mM','Wildtype 400mM'};
for gi = 1:numel(glist)
    gn = glist{gi};
% get variables -----------------------------
x = repmat(timecolname,size(Out.(gn).mean,1),1)';
y = Out.(gn).mean';
e = Out.(gn).se';
% ------------------------------------------


% figure ---------------------------------------
% Create figure
figure1 = figure('Visible','on','Position',[0 0 3 4]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple error bars using matrix input to errorbar
errorbar1 = errorbar(x,y,e,'Marker','o');
dispname = {'Acceleration','Reversal','Pauses'};
color = [[1 0 0]; [0 0.447058826684952 0.74117648601532]; [0 0 0]];
for ei = 1:size(x,2)
set(errorbar1(ei),'DisplayName',dispname{ei},...
        'MarkerFaceColor',color(ei,:),...
    'MarkerEdgeColor',color(ei,:),...
    'Color',color(ei,:));

end

% Create xlabel
xlabel('time(s)');
% Create title
title(titlename{gi});
% Create ylabel
ylabel('P (responses)');
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[95 400]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 1]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','north','EdgeColor',[1 1 1]);

% save
savename = gn;
printfig(savename,pM,'w',5,'h',4)
end
%% ------------------------------------------------

%% plot errorbar across time (manual)
titlename = {'Wildtype 0mM','Wildtype 400mM'};
for gi = 1:numel(glist)
    gn = glist{gi};
% get variables -----------------------------
x = repmat(timecolname,size(Out.(gn).mean,1)-1,1)';
y = Out.(gn).mean(1:2,:)';
e = Out.(gn).se(1:2,:)';
% ------------------------------------------


% figure ---------------------------------------
% Create figure
figure1 = figure('Visible','on','Position',[0 0 3 4]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple error bars using matrix input to errorbar
errorbar1 = errorbar(x,y,e,'Marker','o');
dispname = {'Acceleration','Reversal','Pauses'};
color = [[1 0 0]; [0 0.447058826684952 0.74117648601532]; [0 0 0]];
for ei = 1:size(x,2)
set(errorbar1(ei),'DisplayName',dispname{ei},...
        'MarkerFaceColor',color(ei,:),...
    'MarkerEdgeColor',color(ei,:),...
    'Color',color(ei,:));

end

% Create xlabel
xlabel('time(s)');
% Create title
title(titlename{gi});
% Create ylabel
ylabel('P (responses)');
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[95 400]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 1]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','north','EdgeColor',[1 1 1]);

% save
savename = [gn,' no pauses'];
printfig(savename,pM,'w',5,'h',4)
end
%% ------------------------------------------------
return











% save figure ----------------
pSaveFig = fullfile(pM,'Figure',tname);
if ~isdir(pSaveFig); mkdir(pSaveFig); end
savename = fullfile(pSaveFig, strain );
printfig(savename,pM,'w',3,'h',3)
% -------------------------------------------------

% statistics -------------------
T = table;
T.gname = MWTDBv.groupname;
t = RespNv;
T = [T array2table(t)];
% rid of nan output
T(any(isnan(RespNv),2),:) = [];
[ST,textfile] = RType_stats(T);
% ------------------------------

% save result % ------------------
pSaveFig = fullfile(pM,'Data',tname);
if ~isdir(pSaveFig); mkdir(pSaveFig); end
savename = fullfile(pSaveFig, strain);
gname = sortN2first(gnlist,gnlist);
save(savename, 'R','leg','plateinfo','ST','textfile','Mean', 'SEM', 'SampleSize','gname');
% ---------------------------------

% save text result ---------------
pSaveFig = fullfile(pM,'Stats',tname);
if ~isdir(pSaveFig); mkdir(pSaveFig); end
savename = fullfile(pSaveFig, [strain,'.txt'] );
fid = fopen(savename,'w');
fprintf(fid,'%s',textfile);
fclose(fid);
% -------------------------------   




% report done --------------
fprintf('\n--- Done ---\n\n' );
return
% --------------------------






























