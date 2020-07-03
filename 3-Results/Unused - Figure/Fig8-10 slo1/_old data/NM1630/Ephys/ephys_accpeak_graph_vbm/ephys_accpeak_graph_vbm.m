% function ephys_accpeak_graph2_stats(strainT)

%% INITIALIZING
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
addpath('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/Graphs/ephys');
pM = setup_std(mfilename('fullpath'),'RL','genSave',true);
addpath(pM);
%% paths
pSave = pM;



%% SETTING
pSF = '/Users/connylin/Dropbox/RL Pub PhD Dissertation/Chapters/4-EARS/Data/10sIS by strains';

% Settings
Setting = struct;
Setting.Color = {[0 0 0]; ...
    [1 0 0];...
    [0 0.447058826684952 0.74117648601532];...
    [0.494117647409439 0.184313729405403 0.556862771511078]};
Setting.Marker = {'.' '.' '.' '.'};
Setting.MarkerEdgeColor = Setting.Color;
Setting.MarkerFaceColor = Setting.Color;
Setting.MarkerSize = {2 2 2 2};

rmfactorName = 't';
factor1Name = 'groupname';
pvlimit = 0.001;
alphanum = 0.05;
timeNameList = {'t28_30'; 't1'};

    
%% strains
strainlist = dircontent(pSF);
% get strain info
load('/Users/connylin/Dropbox/Code/Matlab/Library RL/Modules/MWTDatabase/StrainNames.mat');
strainNames(~ismember(strainNames.strain,strainlist),:) = [];


%% load data
for si =1:size(strainNames,1)
    %% get strain info
    strain = strainNames.strain{si}; 
    fprintf('%d/%d: %s\n',si, numel(strainlist), strain);
    genotype = strainNames.genotype{si};
    
    %% graph space out graphs
    statString = '';
    StatOut = struct;
    
    pS = fullfile(pSave,'Stats txt');
    if ~isdir(pS); mkdir(pS); end
    cd(pS);
    fid1 = fopen([strain,' Stats.txt'],'w');

    for ti = 1:numel(timeNameList)
        %% reporting
        timeName = timeNameList{ti};
        fprintf(fid1,'*** %s ***\n',timeName);
        load(sprintf('%s/%s/ephys graph/data_ephys_%s.mat',pSF,strain,timeName));
        
        %% prep graph data
        % get data 
        G = struct;
        % get group names list
        gnlist = {'N2','N2_400mM',strain, [strain,'_400mM']};
        % get group names
        gnu = cell(size(DataG,2),1);
        for gi = 1:size(DataG,2)
           gnu{gi} = DataG(gi).name;
        end
        [i,j] = ismember(gnu,gnlist);
        gnuindseq = j';

        %% organize data
        for gi = 1:numel(gnuindseq)
            gii = gnuindseq(gi);
            x = DataG(gii).time(1,6:41);
            Y = DataG(gii).speedbm(:,6:41);
            Y(:,6) = NaN;
            y = nanmean(Y);
            n = sum(~isnan(Y));
            n = sqrt(n-1);
            n(n<=0) = NaN;
            e = nanstd(Y)./n;
            e2 =e*2;
            g = {DataG(gii).name};
            G.g(1,gi) = g;
            G.x(:,gi) = x';
            G.y(:,gi) = y';
            G.e(:,gi) = e2';
        end
        Setting.DisplayName = G.g;
        
        %% prepare output
        n = size(G.g');
        colNames = celltakeout({strjoinrows([G.g', repmat({'x'},n)], '_');
                    strjoinrows([G.g', repmat({'y'},n)], '_');
                    strjoinrows([G.g', repmat({'e'},n)], '_')})';
        a = [G.x G.y G.e];
        T = array2table(a,'VariableNames',colNames);
        
        pS = fullfile(pSave,['Graph csv ' timeName]);
        if ~isdir(pS); mkdir(pS); end
        cd(pS);
        savename = sprintf('%s graph values.csv',strain);
        writetable(T,savename);
        
        %% graph
        pSavePDF = fullfile(pSave,['Graph ' timeName]);
        if ~isdir(pSavePDF); mkdir(pSavePDF); end
        
        graphEy4(G,genotype,strain,timeName,pSavePDF,Setting,'zeroline',true);

        
        %% STATS
        % get data 
        D = cell(size(DataG,2),1);
        GN = cell(size(DataG,2),1);
        % get time
        time1 = -0.5;
        time2 = 1;
        for gi = 1:size(DataG,2)
            t = DataG(gi).time(1,:);
            ti1 = find(t==time1);
            ti2 = find(t==time2);
            x = DataG(gi).time(1,ti1:ti2);
            Y = DataG(gi).speedbm(:,ti1:ti2);
            Y(:,x==0) = NaN;
%             x(:,x==0) = [];
            D{gi} = Y;
            GN{gi} = repmat({DataG(gi).name},size(Y,1),1);
        end
           
        timeX = x;
        D = cell2mat(D);
        GN= celltakeout(GN);
        
        %% prep for rmanova
        % group name
        a = regexpcellout(GN,'_','split');
        b = a(:,2);
        b(cellfun(@isempty,b)) = {'0mM'};
        s = a(:,1);
        TM = table;
        TM.(factor1Name) = GN;
        TM.strain = s;
        TM.dose = b;
        
        % velocity
        c1 = find(timeX==0)+1;
        c3 = find(timeX==0)-1;
        c2 = find(timeX==time2);
        D1 = D(:,[c3 c1:c2]);
        n = size(D1,2);
        tpu = (0:n-1)';
        a = num2cellstr(tpu);
        t_last = a{end};
        colNames = strjoinrows([cellfunexpr(a,rmfactorName),a],'');
        tptable = table(tpu,'VariableNames',{rmfactorName});

        T = array2table(D1,'VariableNames',colNames);
        
        % anova table
        T = [TM T];
        
        
        
        %% anova
        factorName = 'strain*dose';
        astring = sprintf('%s0-%s%s~%s',rmfactorName,rmfactorName,t_last,factorName);
        rm = fitrm(T,astring,'WithinDesign',tptable);
        StatOut.(timeName).fitrm = rm;
        t = ranova(rm); 
        StatOut.(timeName).ranova = t;
        t = anovan_textresult(t,0, 'pvlimit',pvlimit);
        fprintf(fid1,'RMANOVA(%s,t:%s):\n%s\n',timeName,factorName,t);
       
        % run pairwise
        astring = sprintf('%s0-%s%s~%s',rmfactorName,rmfactorName,t_last,factor1Name);
        rm = fitrm(T,astring,'WithinDesign',tptable);
        StatOut.(timeName).fitrmg = rm;
        t = multcompare(rm,factor1Name);
        StatOut.(timeName).mcomp_g = t;
        t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
        fprintf(fid1,'\nPosthoc(Tukey)curve by group:\n%s\n',t);
        
       
        % comparison by taps
        t = multcompare(rm,factor1Name,'By',rmfactorName);
        StatOut.(timeName).mcomp_g_t = t;
        t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
        fprintf(fid1,'\nPosthoc(Tukey)%s by %s:\n%s\n',factor1Name,rmfactorName,t);

        % comparison within group bewteen time
        t = multcompare(rm,rmfactorName ,'By',factor1Name);
        StatOut.(timeName).mcomp_t_g = t;
        t = multcompare_text2016b(t,'pvlimit',pvlimit,'alpha',alphanum);
        fprintf(fid1,'\nPosthoc(Tukey)%s by %s:\n%s\n',rmfactorName,factor1Name,t);
        fprintf(fid1,'\n\n');
    
        
    end
    fclose(fid1);
    % get n
    StatOut.n_sample = nan(size(DataG,2),1);
    for gi = 1:size(DataG,2)
        StatOut.n_sample(gi) = size(DataG(gi).speedb,1);
    end
    % save
    pS = fullfile(pSave,'mat files');
    if ~isdir(pS); mkdir(pS); end
    cd(pS);
    save(strain,'StatOut');

end


fprintf('Done\n');
% end


























    
    
    
    
    
    
    
    
    
    
    
    