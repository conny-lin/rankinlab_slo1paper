%% INITIALIZING
clc; clear; close all;
%% PATHS
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pM = setup_std(mfilename('fullpath'),'RL','genSave',false);
pSave = fileparts(pM);

%% get data
% load([fileparts(pSave),'/Dance_DrunkMoves/Dance_DrunkMoves.mat'])
msrlist = {'RevDur','RevFreq','RevSpeed'};
gnu = {'N2','N2_400mM'};
isiu = {'10','60'};
posthoctype = 'bonferroni';
T1 = table;

cd(pSave); 
fid1 = fopen(sprintf('MANOVA %s.txt',posthoctype),'w');
for msri = 1:numel(msrlist)
    D = cell(numel(isiu)*numel(gnu),1);
    G = D;
    G2 = D;
    n = 1;
    for isii = 1:numel(isiu)
    for gni = 1:numel(gnu)
        isi = isiu{isii};
        pD = sprintf('%s/Data/Hab curve %ssISI recovery/Dance_Glee_Acafellas/Data Raw',...
            fileparts(pSave),isi);
        msr = msrlist{msri};
        gn = gnu{gni};
        clear p;
        p = sprintf('%s/%s_Mean %s.csv',pD,msr,gn);
        T = readtable(p);
        % calculate % recovery
        d = T.s31./T.s1;
        g = repmat(cellstr([{gn},{isi}]),numel(d),1);
        g2 = strjoinrows(g);
 
        
        % add to master
        D{n} = d;
        G{n} = g;
        G2{n} = g2;
        n = n+1;
    end
    end
    
    % take out cell
    G = celltakeout(G,'multirow');
    G2 = celltakeout(G2,'multirow');
    G1 = cell(1,size(G,2));
    for c = 1:size(G,2)
       G1{c} = G(:,c);
    end
    D = cell2mat(D);
    G = G1; 
    
    %% summary table
    T = grpstatsTable(D,G2);
    T = [cell2table(repmat({msr},size(T,1),1),'VariableName',{'msr'}) T];
    T1 = [T1; T];
    
    %% anova
    fprintf(fid1,'-- %s --\n',msr);
    close all;
    gvarnames = {'isi','dose'};
    [~,t] = anovan(D,G,'varnames',gvarnames,'model','full','display','off');
    result = anovan_textresult(t,'off');
    fprintf(fid1,'MANOVA\n');
    fprintf_cell(result,fid1);
    
    %% posthoc
    
    
    [~,~,s] = anova1(D,G2,'off');
    [c,~,~,gn]= multcompare(s,'display','off','Ctype',posthoctype);
    result = multcompare_text2016b(c,'grpnames',gn);
    fprintf(fid1,'posthoc(%s)\n',posthoctype);
    fprintf(fid1,'%s\n',result);
    
end
clear G1 c n T p gn msr isi pD gni msri d g gnu isii isiu;


fclose(fid1);
   
 writetable(T1,'graph.csv');
   
   
   
   
   
   




