%% calculate stats for recovery
% used posthoc tukey-kramer as conservative assay and accounts for uneven
% sample sizes
%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GLOBAL VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pHome = fileparts(fileparts(mfilename('fullpath')));
pHome = '/Users/connylin/Dropbox/Publication/Manuscript RL Alcohol hab model/Figures Tables Data/Fig3 Spontaneous recovery';
ctype = 'tukey-kramer';
msrlist = {'RevDur','RevFreq','RevSpeed'};
grouplist = {'N2','N2_400mM'};
ISI = [10 60];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% calculate recovery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fullfile(pSave,'MANOVA.txt'),'w'); % open text output file
for msri = 1:numel(msrlist) % cycle by measures
    % create output array
    L = {};
    D = [];
    for x = 1:numel(ISI) % cycle through isi
        for gi = 1:numel(grouplist) % cycle through groups
            msr = msrlist{msri}; % get msr name
            gn = grouplist{gi}; % get group name
            isi = ISI(x); % get isi name
            foldername = sprintf('%dsISI recovery',isi);
            str = fullfile(pHome, foldername, 'Dance_Glee_Acafellas','Data Raw',sprintf('%s_Mean %s.csv',gn)); % get path to raw data
            a = readtable(str); % read table
            % calculate hab level
            hab = mean([a.s28 a.s29 a.s30],2);
            initial = a.s1;
            rec = a.s31;
            % difference between hab to initial
            habp = hab./initial;
            recp = (rec./initial)- habp;
            % cal stats
            leg = [repmat({num2str(isi)},numel(recp),1) repmat({gn},numel(recp),1)];
            % add to main array 
            L = [L;leg];
            D = [D;recp];
        end
    end
    
    % create table
    T = table;
    T.ISI = L(:,1);
    T.group = L(:,2);
    T.y = D;
    
    % anovan
    return
    anovan_std(Y,G,gvarnames,pSave,prefix,varargin)

    
    % anovan
    y = T.y;
    group = {T.ISI T.group};    
    [p,t,s,tr] = anovan(y,group,'varnames',{'ISI','dose'},...
                        'display','off','nested',[0 0 ; 1 0]);
    str = anovan_textresult(t);
    fprintf(fid,'%s -----------------------\n',msr);
    fprintf(fid,'MANOVA ISI x dose (nested)\n');
    for x = 1:numel(str)
        fprintf(fid,'%s\n',str{x});
    end
    
    %% t test of significant recovery
    % The alternative hypothesis is that the population distribution does 
    % not have a mean equal to zero. 
    fprintf(fid,'\nt test for recovery > zero\n');
    isiu = unique(T.ISI);
    gu = unique(T.group);
    for isii = 1:numel(isiu)
        for gi = 1:numel(gu)
            isi = isiu{isii};
            gn = gu{gi};
            i = ismember(T.ISI,isi) & ismember(T.group,gn);
            d = T.y(i);
            [h,p,ci,s] = ttest(d);
            if p < 0.0001
                pstr = sprintf('< 0.0001');
            else
                pstr = sprintf('= %.4f',p);
            end
            fprintf(fid,'%sISI(%s): %s\n',isi,gn,pstr);
        end
    end
    
    %% posthoc between dose within ISI
    fprintf(fid,'\nposthoc %s (dose)\n',ctype);
    a = unique(T.ISI);
    for ai = 1:numel(a)
        i = ismember(T.ISI,a{ai});
        d = T.y(i);
        g = T.group(i);
        [~,~,s] = anova1(d,g,'off');
        [c,~,~,gnames] = multcompare(s,'display','off','ctype',ctype);
        [~,str] = multcompare_pairinterpretation(c,gnames);
        for x = 1:numel(str);
            fprintf(fid,'%sISI: %s\n',a{ai},str{x});
        end
    end
    
    %% posthoc between ISI within dose
    fprintf(fid,'\nposthoc %s (ISI)\n',ctype);
    a = unique(T.group);
    for ai = 1:numel(a)
        i = ismember(T.group,a{ai});
        d = T.y(i);
        g = T.ISI(i);
        for x = 1:numel(g)
           g{x} = sprintf('ISI%s',g{x}); 
        end
        [~,~,s] = anova1(d,g,'off');
        [c,~,~,gnames] = multcompare(s,'display','off','ctype',ctype);
        [~,str] = multcompare_pairinterpretation(c,gnames);
        for x = 1:numel(str);
            fprintf(fid,'%s: %s\n',a{ai},str{x});
        end
    end
    
    %% next line
    fprintf(fid,'\n\n');

end
fclose(fid); close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


