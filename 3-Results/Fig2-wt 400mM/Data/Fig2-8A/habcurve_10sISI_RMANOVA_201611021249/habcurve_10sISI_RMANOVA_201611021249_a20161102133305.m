%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% GLOBAL INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths ++++++++++++++
pD = fileparts(pM);
pSave = pM;
% --------------------


% settings ++++++
% ---------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data ================================================================
% IMPORT HABITUATION CURVES DATA +++++++++++
% load MWTDB 
load(fullfile(fileparts(pD),'MWTDB.mat'),'MWTDB');
pMWT = MWTDB.mwtpath;
% get hab data
[Data,msrlist] = importHabCurveData(pMWT);
% ------------------------------------------

% create tap data obj ++++++++++++++++++++++
TD = TapData;
v = {'mwtpath','groupname','tap'};
for vi = 1:numel(v)
    TD.(v{vi}) = Data.(v{vi});
end
TD.Response = table2array(Data(:,msrlist));
TD.Measures = msrlist;
% ------------------------------------------
% =========================================================================



% anova ===================================================================
rmanovar = rmANOVA(TD);
for msri = 1:numel(msrlist)
    
    msr = msrlist{msri};
    textout = rmanovar.(msr);
    % save text data +++++++
    p = fullfile(pSave,sprintf('%s RMANOVA.txt',msr));
    fid = fopen(p,'w');
    fprintf(fid,'%s',textout);
    fclose(fid);
    % ----------------------

end
% =========================================================================


% habituated level ========================================================
D = pctData(TD,1,30);
% remove invalid data
D((isnan(D.pct) | isinf(D.pct)),:) = [];

% generate stats
D.gname = regexprep(D.gname,'N2_400mM','400mM');
D.gname = regexprep(D.gname,'N2','0mM');
x = D.pct;
g = D.gname;
gnu = unique(g);

[text,T,p,s,t,ST] = anova1_autoresults(x,g);
T.Properties.RowNames = T.gnames;
str = sprintf('*** Descriptive stats ***');
gn = strjoin(gnu',', ');
str = sprintf('%s\ngroupnames = %s',str,gn);
cmp = {'N','mean','SE'};
for ci = 1:numel(cmp)
    c = cmp{ci};
    n = num2cellstr(T.(c)');
    n2  = strjoin(n,', ');
    str = sprintf('%s\n%s = %s',str,c,n2);
end

[txt] = anova1_multcomp_auto(x,g);
str = sprintf('%s\n\n*** ANOVA ***\n%s',str,txt);


y = 1;
str = sprintf('%s\n\n*** t-test against %d ***',str,y);
for gi = 1:numel(gnu)
    gg = gnu{gi};
    i = ismember(g,gg);
    xx = x(i);
    text = ttest_auto(xx,y);
    
    
    str = sprintf('%s\n%s, %s',str,gg,text);
end

% save text data +++++++
p = fullfile(pSave,sprintf('%s hab level ttest.txt',msr));
fid = fopen(p,'w');
fprintf(fid,'%s',str);
fclose(fid);
% ----------------------
% =========================================================================
% 
% 
% 
% 
% %% Graph ==================================================================
% % plot line graph +++++++++++++++++++++++++++
% w = 4;
% h = 3;
% for msri = 1:numel(msrlist)
%     msr = msrlist{msri};
%     
%     fig1 = plotHabStd(TD,msr,'off');
%     set(fig1,'Visible','on')
% 
%     % adjust color +++++
%     color = color_rainbow;
%     f1 = get(fig1);
%     a1 = get(f1.Children);
%     e1 = get(a1.Children);
%     % apply settings 
%     gp = graphsetpack('bkrainbowline');
%     gpn = fieldnames(gp);
%     for gi = 1:numel(gp.Color)
%         for gpi = 1:numel(gpn)
%             fig1.Children.Children(gi).(gpn{gpi}) = gp.(gpn{gpi}){gi};
%         end
%     end
% 
%     % ------------------
%     
%     % save +++++++
%     savename = sprintf('%s/%s',pSave,msr);
%     printfig(savename,pSave,'w',w,'h',h,'closefig',1);
%     % ------------
% end
% % ------------------------------------------
% % =========================================================================


%% Save and exit ==========================================================
% save data ++++++++++++
save(fullfile(pSave,'data.mat'),'Data','TD','rmanovar');
% ----------------------

% report done ++++++++++++
fprintf('\nDONE\n'); return
% ------------------------
% =========================================================================

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















