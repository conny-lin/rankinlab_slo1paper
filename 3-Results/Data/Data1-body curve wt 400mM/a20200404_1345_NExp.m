%% analyze body curve for all strans.
% 2020-04-04 13:46
%% paths (program)
clc
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General')
add_code_package_paths('RL')

%% Define save paths
pSave = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Data1-body curve wt 400mM';

%% get experiment information 
pData = '/Users/connylin/Dropbox/CA/CA Publications/Manuscript RL Alcohol hab model slo1/3-Results/Data1-body curve wt 400mM/Dance_InitialEtohSensitivityPct/data.mat';
load(pData);

%% calculate stats by exp
expname = unique(MWTDB.expname);
gname = unique(MWTDB.groupname);
rown = numel(expname);
coln = numel(gname);

MeanByExp = nan(rown,coln);

for ei = 1:numel(expname)
    e = expname(ei);
    for gi = 1:numel(gname)
        g = gname(gi);
        i = ismember(MWTDB.expname, e) & ismember(MWTDB.groupname, g);
        mwtid = MWTDB.mwtid(i);
        d = DataMeta.curve(ismember(DataMeta.mwtid,mwtid));
        dmean = nanmean(d);
        MeanByExp(ei,gi) = dmean;  
    end
end

%% ETL
data = reshape(MeanByExp,numel(MeanByExp),1);
g = repmat(gname',rown,1);
group = reshape(g,numel(g),1);
%% run stats
[txt,anovastats,multstats,T,ST,ST2] = anova1_std_v2(data,group);
%% save stats
filesavename = fullfile(pSave,'ANOVA curve.txt');
fid = fopen(filesavename,'w');
fprintf(fid,'%s',txt); 
fclose(fid);
%% save table
pS = fullfile(pSave,'curve.csv');
writetable(T,pS)




%% analysis
% Dance_InitialEtohSensitivityPct(pMWT,'pSave',pSave);
% Dance_Raster2(pMWT,'pSave',pSave);

