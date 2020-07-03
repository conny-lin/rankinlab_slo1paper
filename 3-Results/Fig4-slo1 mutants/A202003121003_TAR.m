%% analyze rType for all strans.
% 2020-02-16 12:15
%% add paths (program)
clear
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% set source and save paths
% p = mfilename('fullpath');
% [~,f] = fileparts(p);
% pSave = create_savefolder(p, f);
pSave = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Figures & Data/Fig4 slo1 mutants';


%% define by group
StrainNameTarget = {'JPS429','BZ142','NM1968','NM1630'};


%% TAR

fprintf('TAR*****\n');

for i = 1:numel(StrainNameTarget)
    % get strain name
    strainname = StrainNameTarget{i};
    % display running
    fprintf('analyzing %s\n',strainname);
    
    % load data
    p = fullfile(pSave,strainname,'TAR','Dance_rType','data.mat');
    load(p)
    fprintf('got data\n');
    
    
    %% transform Dance_rType into excel output
    % restructure data
    d = [G.AccProb.X(:,1)';...
        G.AccProb.Y'; ...
        G.AccProb.E'; ...
        G.AccProb.N'];
    % make label name
    g = repmat(G.AccProb.groupname,3,1);
    g = [{'tap'};g];
    % create stats name
    gn = numel(G.AccProb.groupname);
    stats = [{'tap'};...
            repmat({'mean'},gn,1);...
            repmat({'se'},gn,1);...
            repmat({'n'},gn,1)];
    %% create table column names
    t = G.AccProb.X(:,1);
    t = num2cell(t);
    t = cellfun(@num2str,t,'UniformOutput',0);
    tname = repmat({'t'},numel(t),1);
    for x = 1:numel(t)
       a = strjoin([tname(x),t(x)],'');
       a = cellstr(a);
       tname(x) = a;
    end

    %% create table
    T = table;
    T.groupname = g;
    T.statsname = stats;
    for x = 1:numel(tname)
       T.(tname{x})= d(:,x); 
    end
    % write table
    fname = fullfile(pSave,strainname,'TAR','Dance_rType','AccProb_patch.csv');
    writetable(T,fname)

end










