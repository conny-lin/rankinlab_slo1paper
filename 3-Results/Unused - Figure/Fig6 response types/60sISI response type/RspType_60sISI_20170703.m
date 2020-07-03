%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GLOBAL INFORMATION
% settings %--------------------
n_lowest = 10;
tapstarttime = 100;
isi = 60;
tapN = 30;
assaytimeperiod = 1; % 1s after tap
baselinetimeperiod = 0.3; % 0.3 seconds before tap
responseTypes = {'accelerate forward'; 'reversal';'pause'};
rTarget = {'accelerate forward'; 'reversal'};
% construct times
taptimelist = [tapstarttime:isi:(tapstarttime+isi*(tapN-1))];
matfiletimelist = taptimelist-1;
% --------------------------------


%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pData = fullfile(pSave,'Dance_Raster2');
if ~isdir(pData)
    Data_Raster = Dance_Raster2(pMWT,'pSave',pSave_Raster,'graphopt',false); % Dance raster    
end
% strains ------
[~,~, glist] = dircontent(pData);
glist(regexpcellout(glist,'graffle')) = [];
%----------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        if responsei+5 > assaytimei(end)
            r2 = assaytimei(end);
        else
            r2 = responsei+5;
        end
        responsetimei = responsei:1:r2;
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
    

    
    % organize descripive output data
    Out.(groupname).n = NVAL;
    Out.(groupname).mean = MEAN;
    Out.(groupname).se = SEM;
    
    % table output for descriptive data
    T1 = statsTable_multvar(NVAL,MEAN,SEM,responseTypes);
    T = table;
    T.tap = (1:numel(taptimelist))';
    T = [T T1];
    writetable(T,fullfile(pSave,sprintf('%s.csv',glist{gi})));
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
save(fullfile(pSave,'data.mat'),'Out','DataStats','responseTypes');



%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings
rmName = 'tap';
factors = {'group','response_type'};
dvname = 'probability';
idname = 'mwtid';
alpha = 0.05;
pvlimit = 0.001;

% prep data
D = DataStats;
D(~ismember(D.response_type,rTarget),:) = []; % delete irrelevant data
% simplify response type name
D.response_type = regexprep(D.response_type,'accelerate forward','acc');
D.response_type = regexprep(D.response_type,'reversal','rev');
[textout,DS] = anovarm_std(D,rmName,factors,idname,dvname); % rmANOVA
% export to text
fid = fopen(fullfile(pSave,'RMANOVA.txt'),'w');
fprintf(fid,'%s',textout);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    fig_pResp(pSave,x,y,e,titlename{gi},gn,'dispname',rTarget); % save figure
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% report done --------------
fprintf('\n--- Done ---\n\n' );
return

%% --------------------------








%%
A = DS.posthoc.factors_by_tap
i = A.pValue < 0.05 & ismember(A.factors_1,'N2_400mM*acc') & ismember(A.factors_2,'N2_400mM*rev')
a = A(i,{'factors_1','factors_2','tap','pValue'});
b = a.tap;
strjoin(num2cellstr(b)',', ')


b = a.pValue;

strjoin(num2cellstr(b)',', ')





















