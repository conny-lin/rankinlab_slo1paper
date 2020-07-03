%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SETTING
genotype = 'wild-type'; % get genotype
figname = 'N2';
w = 9;
h = 4;

%% GET MWT INFORMATION
% get mwt paths
pDatabase = fullfile('/Volumes/COBOLT/MWT','MWTDB.mat');
load(pDatabase);
MWTDBM = MWTDB.text;
% select data
i = ismember(MWTDBM.rc,'100s30x60s10s')...
    & ismember(MWTDBM.groupname,{'N2_400mM'})...
    & MWTDBM.exp_date > 20120301 ...
    ;
e = MWTDBM.expname(i);
i = ismember(MWTDBM.expname,e) & ismember(MWTDBM.groupname,{'N2','N2_400mM'});

MWTDB = MWTDBM(i,:);
pMWT = MWTDB.mwtpath;
fprintf('\nexp report: \n');
tabulate(MWTDB.groupname)
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% DANCE ANALYSIS COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL SENSITIVITY 
fprintf('\n\n --- CURVE ANALYSIS ---\n');
Data_Sen = Dance_InitialEtohSensitivityPct_v1707(pMWT,pSave); % run initial sensitivity
pData_Sen = Data_Sen.pSave;
A = Data_Sen.Stats.dscstat;
Summary.Curve.groupname = gn;
Summary.Curve.N = A.N;
Summary.Curve.Y = A.mean;
Summary.Curve.E = A.SE;

fprintf('\n\n --- CURVE ANALYSIS DONE ---\n');
% --------------------------


% REVERSAL
fprintf('\n\n --- REVERSAL ANALYSIS ---\n');
Data_SS = Dance_ShaneSpark5(pMWT,'pSave',pSave);
pData_SS = Data_SS.PATHS.pSaveA;
load(fullfile(pData_SS,'Dance_ShaneSpark5.mat')); % load data
D = MWTSet.ByGroupPerPlate.RevFreq;
R = TWR_graph_data(D);
Summary.TWR = R;
fprintf('\n\n --- REVERSAL ANALYSIS DONE ---\n');
% -------------------------------


% ACCELERATION
fprintf('\n\n --- ACCELERATION ANALYSIS ---\n');
Data_TAR = Dance_rType2(pMWT,'pSave',pSave,'genotype',genotype); % rType
pData_TAR = Data_TAR.PATHS.pSaveA;
A = load(fullfile(pData_TAR,'data.mat'));
Summary.TAR = A.G.AccProb;
fprintf('\n\n --- ACCELERATION ANALYSIS DONE ---\n');

% -----------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markersize = 4.5;
gp = graphsetpack('cathyline');
gpn = fieldnames(gp);

close all
figure1 = figure('Visible','off');

% bar graph +++++++++++++
GN = Summary.Curve.groupname;
GN = regexprep(GN,' 400mM','');
if ~strcmp(GN{1},'N2'); error('bad n2');end
Y = Summary.Curve.Y;
E = Summary.Curve.E;
% Create axes
axes1 = axes('Parent',figure1,'Position',[0.07 0.11 0.1 0.45],'Box','off');
hold(axes1,'on');
% axes properties
set(axes1,'XTick',[1 2],'XTickLabel',{'+','-'});
xlim([0.5 2.5]);
ylim([-60 20]);
ylabel('% curve (etoh/control)')
bar1 = bar(Y,'Parent',axes1,'EdgeColor',[0 0.447058826684952 0.74117648601532],'FaceColor',[0 0.447058826684952 0.74117648601532]);
% errorbar
e1 = errorbar(Y,E,'LineStyle','none','Color','k','LineWidth',1,'Marker','none');
% ------------------------


% plot 2 ++++++++++++++++
GN = Summary.TWR.gn;
X = Summary.TWR.X;
Y = Summary.TWR.Y;
E = Summary.TWR.E;
if ~strcmp(GN{1},'N2')
    X = sortN2first(GN,X')';
    Y = sortN2first(GN,Y')';
    E = sortN2first(GN,E')';
    GN = sortN2first(GN,GN);
end 
% Create axes
axes2 = axes('Parent',figure1,'Position',[0.25 0.11 0.3 0.8]);
hold(axes2,'on');
xlim([0 30.5]);
ylim([0,1]);
ylabel('P (reversal)');
xlabel('Tap');
e2 = errorbar(X,Y,E,'Marker','o','MarkerSize',markersize);
e2 = graphApplySetting(e2,'cathyline');
% Create multiple error bars using matrix input to errorbar
legname = {'+ 0mM','+ 400mM','(-) 0mM','(-) 400mM'};
for ei = 1:numel(GN)
    set(e2(ei),'DisplayName',legname{ei});
end
% Create legend
legend1 = legend(axes2,'show');
set(legend1,'Position',[0.04 0.66 0.1 0.22],'EdgeColor',[1 1 1]);
annotation(figure1,'textbox',[0.02 0.95 0.08 0.08],...
    'String',{genotype},'FontWeight','bold','EdgeColor','none');
% ------------------------


% plot 3 +++++++++++
GN = Summary.TAR.groupname;
X = Summary.TAR.X;
Y = Summary.TAR.Y;
E = Summary.TAR.E;
if ~strcmp(GN{1},'N2')
    X = sortN2first(GN,X')';
    Y = sortN2first(GN,Y')';
    E = sortN2first(GN,E')';
    GN = sortN2first(GN,GN);
end    
% Create axes
axes3 = axes('Parent',figure1,'Position',[0.65 0.11 0.3 0.8]);
hold(axes3,'on');
xlim([0 30.5])
 ylim([0,1]);

ylabel('P (acceleration)');
xlabel('Tap');
e3 = errorbar(X,Y,E,'Marker','o','MarkerSize',markersize);
e3 = graphApplySetting(e3,'cathyline');
% ------------------------

% save -------------------------------------------
printfig(figname,pSave,'w',w,'h',h,'closefig',1);
% -------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% STATS REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% report_All(pData_Sen,pData_TAR,pData_SS,pSave,MWTDB,figname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n ~~~~ DONE ~~~~\n');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
