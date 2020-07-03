
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
timelist = [95:10:29*10+95];
time_baseline = [-0.3 -0.1];
% time_response = [0.1 0.5];
% statsName = 'exclude no response';
% 
% % create time frames
% startList = [98 388];
% endList = startList+assaywindow;
% --------------------------------


%% RESPONSE PROBABILITY: ACC, PAUSE, NO RESPONSE, REVERSAL, DECELLEARTION
% construct output array -------------
a = 0:0.2:1;
b = repmat(a',1,numel(timelist));
c = repmat((timelist+5),numel(a),1);
d = b+c;
timecolname = reshape(d,numel(d),1)';
ncol = numel(timecolname);
nrow = numel(glist);
Out = struct;
% -----------------------------

    
for gi = 1:numel(glist) % cycle through groups
    
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
    responseTypes = {'accelerate forward'; 'reversal';'pause'};
    MEAN = nan(numel(responseTypes), ncol);
    SEM = MEAN;
    NVAL = MEAN;
    % -----------------------------------
    
    for ti = 1:numel(timelist) % cycle through time choices

        % report progress %------------
        fprintf(''); % separator
        processIntervalReporter(numel(timelist),1,'time',ti);
        % ------------------------------

        % time info -----------------
        time1now = timelist(ti);
        taptimenow = time1now + 5;
%         tname = timelist{ti};
%         tfilename = ['data_ephys_',tname,'.mat'];
        % --------------------------

        % data % -------------
        f = matfileinfo.fname{matfileinfo.t1 == time1now};
        p = fullfile(pData, groupname,f);
        load(p);
        % find tap time column
        i = find(rTime(1,:) >= taptimenow & rTime(1,:) <= taptimenow+1);
        i = i(1);
        % get data before tap Baseline
        Baseline = Data(:,i-3:1:i-1);
        Response = Data(:,i:1:i+5); % 10 s after the tap
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
            
            % transform data 
            
            return
            %% descriptive stats 
            [d,leg] = statsBasic(RespNv);
            m = d(:,3);
            se = d(:,5);
            n = d(:,1);
            t = taptimenow:0.2:(taptimenow+(5*0.2));
            
            
            % place data in output array
            [i,j] = ismember(t,timecolname);
            MEAN(ri,j) = m;
            SEM(ri,j) = se;
            NVAL(ri,j) = n;
            % ------------------------------------
            return
            
        end 
    end
    
    Out.(groupname).mean = MEAN;
    Out.(groupname).se = SEM;
    Out.(groupname).n = NVAL;
end

cd(pM);
save('matlab.mat','Out');



%% FIGURES
% settings %--------------------
% timelist = [95:10:29*10+95];
timelist = 100:10:(100+(29*10));
coli = [1:6:(30*6)];
% --------------------------------

% extract only 0.1s data
glist = fieldnames(Out);
for gi = 1:numel(glist)
    gn = glist{gi};
end
% --------------------------------

% plot errorbar across time (manual)
titlename = {'Wildtype 0mM','Wildtype 400mM'};
for gi = 1:numel(glist)
    gn = glist{gi}; % get group specific info 
    d = Out.(gn); % get data 
    % get variables -----------------------------
    y = d.mean(1:2,coli)';
    e = d.se(1:2,coli)';
    x = repmat(timelist,size(d.mean,1)-1,1)';
    % ------------------------------------------
    fig_pResp(pM,x,y,e,titlename{gi},gn); % save figure
end


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






























