function MWTSet = Dance_ShaneSpark4(pMWT,varargin)
%% INFORMATION



%% DEFAULTS
% MWTSet = struct;

timestamp = 'off';
pSave = '/Users/connylin/Dropbox/RL/Dance Output';
Varin= varargin;
msrlist = {'RevFreq','RevSpeed','RevDur'};
vararginProcessor;


%% standard Dance processing
[MWTSet,pSave] = DanceM_MWTSetStd3(pMWT,mfilename('fullpath'),Varin,pSave,'timestamp',timestamp);

%% get trv data
MWTSet.Raw = import_trv_table(pMWT);


%% make group graph data
DataTrv = MWTSet.Raw;
MWTDB = MWTSet.MWTDB;
[i,j] = ismember(DataTrv.mwtid,MWTDB.mwtid);
DataTrv.groupname = MWTDB.groupname(j(i));
DataTrv.strain = MWTDB.strain(j(i));

gnu = output_sortN2first(unique(DataTrv.groupname)); % group name
tapu = unique(DataTrv.tap); % get unique tap


%% group eStat
G = struct;
header = struct;
header.tap = repmat(tapu,numel(gnu));
for gi = 1:numel(gnu)
for msri = 1:numel(msrlist)
    i = ismember(DataTrv.groupname,gnu(gi));
    d = table2array(DataTrv(i,msrlist{msri}));
    g = table2array(DataTrv(i,{'tap'}));
    T = grpstatsTable(d,g,'gnameutitle','tap');

    G.(msrlist{msri}).groupname(gi) = gnu(gi);
    G.(msrlist{msri}).tap(:,gi) = T.tap;
    G.(msrlist{msri}).Mean(:,gi) = T.mean;
    G.(msrlist{msri}).N(:,gi) = T.n;
    G.(msrlist{msri}).SE(:,gi) = T.se;

end
end
MWTSet.ByGroupPerPlate = G;


%% Graph: ByGroupPerPlate
DataG = MWTSet.ByGroupPerPlate;

% graph mean curve 

for msri = 1:numel(msrlist);
    gn = regexprep(DataG.(msrlist{msri}).groupname,'_',' ');
    y = DataG.(msrlist{msri}).Mean;
    x = DataG.(msrlist{msri}).tap;
    e = DataG.(msrlist{msri}).SE;
    msr = msrlist{msri};
    Graph_HabCurveSS(x,y,e,gn,msr,pSave)
end

%% save mat
cd(pSave);
save([MWTSet.ProgramName,'.mat'],'MWTSet');

return




%% STATISTICS FOR GRAPHS (by group by plate)
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
B = struct;
B.GroupNames = gnames;
for g = 1:numel(gnames)
    gname = gnames{g};

    D1 = D.(gname);
    
    plateN = numel(D1.MWTplateID);
    X = (1:size(D1.N_TotalN,1))';

    B.PlateN(g,1) = plateN;

    B.N_Total.X(:,g) = X;
    B.N_Total.Y(:,g) = nanmean(D1.N_TotalN,2);
    B.N_Total.E(:,g) = nanstd(D1.N_TotalN')';

    B.N_Reversed.X(:,g) = X;
    B.N_Reversed.Y(:,g) = nanmean(D1.N_Reversed,2);
    B.N_Reversed.E(:,g) = nanmean(D1.N_Reversed,2);

    B.N_NoResponse.X(:,g) = X;
    B.N_NoResponse.Y(:,g) = nanmean(D1.N_NoResponse,2);
    B.N_NoResponse.E(:,g) = nanstd(D1.N_NoResponse')';

    B.RevFreq.X(:,g) = X;
    B.RevFreq.Y(:,g) = nanmean(D1.RevFreq_Mean,2);
    B.RevFreq.SD(:,g) = nanstd(D1.RevFreq_Mean')';
    B.RevFreq.E(:,g) = B.RevFreq.SD(:,g)./...
        sqrt(repmat(plateN,size(B.RevFreq.SD(:,g),1),1));

    B.RevDur.X(:,g) = X;
    B.RevDur.Y(:,g) = nanmean(D1.RevDur_Mean,2);
    B.RevDur.SD(:,g) = nanstd(D1.RevDur_Mean')';
    B.RevDur.E(:,g) = B.RevDur.SD(:,g)./...
        sqrt(repmat(plateN,size(B.RevDur.SD(:,g),1),1));

    B.RevSpeed.X(:,g) = X;
    B.RevSpeed.Y(:,g) = nanmean(D1.RevSpeed_Mean,2);
    B.RevSpeed.SD(:,g) = nanstd(D1.RevSpeed_Mean')';
    B.RevSpeed.E(:,g) = B.RevSpeed.SD(:,g)./...
        sqrt(repmat(plateN,size(B.RevSpeed.SD(:,g),1),1));

end
MWTSet.Graph.HabCurve = B;



%% CALCULATION: LAST 3 TAPS 
pSaveA = [pSave,'/','Stats HabLevel'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
M = fieldnames(MWTSet.Graph);
M = {'RevFreq', 'RevDur', 'RevSpeed'};

% calculate data - hab level
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;
for m = 1:numel(M);% for each measure
    B = []; G = {};
    for g = 1:numel(gnames)
        d = D.(gnames{g}).([M{m},'_Mean'])(end-2:end,:);
        % calculate last 3 taps average per plate
        d1 = nanmean(d);
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];
    end
    % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 

    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        % export anova
        vname = regexprep(t(1,:),'Prob>F','P_value');
        T = cell2table(t(2:end,:),'VariableNames',vname);
        cd(pSaveA); 
        writetable(T,[M{m},' HabLevel ANOVA.csv'],'Delimiter',',');
        
        % posthoc
        A.(M{m}).posthoc = multcompare_pairinterpretation(c,gnames);
        % export posthoc
        T = cell2table(a);
        cd(pSaveA); 
        writetable(T,[M{m},' HabLevel posthoc bonferroni.csv'],'Delimiter',',');
        
    end
end
MWTSet.Graph.HabLevel = A;


%% CALCULATION: INITIAL
pSaveA = [pSave,'/','Stats Initial'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;
for m = 1:numel(M);% for each measure
    B = []; G = {};
    for g = 1:numel(gnames)
        d1 = D.(gnames{g}).([M{m},'_Mean'])(1,:);
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];

    end
    
    % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 

    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        % export anova
        vname = regexprep(t(1,:),'Prob>F','P_value');
        T = cell2table(t(2:end,:),'VariableNames',vname);
        cd(pSaveA); 
        writetable(T,[M{m},' Initial ANOVA.csv'],'Delimiter',',');
        
        % posthoc
        A.(M{m}).posthoc = multcompare_pairinterpretation(c,gnames);
        % export posthoc
        T = cell2table(a);
        cd(pSaveA); 
        writetable(T,[M{m},' Initial posthoc bonferroni.csv'],'Delimiter',',');
        
    end
end
MWTSet.Graph.Initial = A;


%% CALCULATION HAB RATE: Andrew's rate - half life
% = first tap achieved half of hab level 
% (hab level (mean of last 3 taps) - initial)/2
M = {'RevFreq', 'RevDur','RevSpeed'};
pSaveA = [pSave,'/','Stats Hab rate slope to half life'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;

for m = 1:numel(M);% for each measure
    % calculate hab rate
    B = []; G = {};
    for g = 1:numel(gnames)
        data = D.(gnames{g}).([M{m},'_Mean']);
        % get initial
        init = data(1,:);
        % get hab level
        hab = nanmean(data(end-2:end,:)); 
        a = data - repmat(hab,size(data,1),1); % minus hab level
        b = a./repmat(a(1,:),size(a,1),1); % divided by initial
        i = b <= .5; % find 50% point
        midtap = nan(size(i,2),1);
        for x = 1:size(i,2)
            a = find(i(:,x));
            if isempty(a) == 0
                midtap(x) =  a(1);
            end
        end
        % get data
        B = [B;midtap];
        n = numel(midtap);
        G = [G; repmat(gnames(g), n,1)];
    end
    
    % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 
    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        % export anova
        vname = regexprep(t(1,:),'Prob>F','P_value');
        T = cell2table(t(2:end,:),'VariableNames',vname);
        cd(pSaveA); 
        writetable(T,[M{m},' Initial ANOVA.csv'],'Delimiter',',');
        
        % posthoc
        A.(M{m}).posthoc = multcompare_pairinterpretation(c,gnames);
        % export posthoc
        T = cell2table(a);
        cd(pSaveA); 
        writetable(T,[M{m},' Initial posthoc bonferroni.csv'],'Delimiter',',');
        
    end
end
MWTSet.Graph.HabRate_halflife = A;



%% CALCULATION HAB RATE: AREA BELOW THE CURVE (INTEGRAL)
% makes no assumption on which is the habituated level
% find the lowest point of response and set it as zero
% calculate area under the curve
% the smaller the area, the faster the animal reaches lowest level
% = first tap achieved half of hab level 
% (hab level (mean of last 3 taps) - initial)/2
M = {'RevFreq', 'RevDur','RevSpeed'};
pSaveA = [pSave,'/','Stats Hab rate integral'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
gnames = fieldnames(D);
A = struct;

for m = 1:numel(M);% for each measure
    % calculate hab rate
    B = []; G = {};
    for g = 1:numel(gnames)
        data = D.(gnames{g}).([M{m},'_Mean']);
        % get lowest response
        datadiff = data - repmat(min(data),size(data,1),1);
        % find area under the curve
        platearea = nan(size(datadiff,2),1);
        for mwti = 1:size(datadiff,2) 
            a = datadiff(:,mwti);
            dint = nan(size(a,2)-1,1);
            for ti = 1:size(a,1)-1
                n1 = a(ti);
                n2 = a(ti+1);
                if n1 > n2; 
                    nmax = n1; 
                    nmin = n2;
                elseif n1==n2; 
                    nmax = n1;
                    nmin = n2;
                else
                    nmax = n2;
                    nmin= n1;
                end
                basearea = nmin*1;
                trianglearea = (nmax*1)/2;
                totalarea = basearea + trianglearea;
                dint(ti) = totalarea;      
            end
            platearea(mwti) = sum(dint);
        end
                
        % get data
        B = [B;platearea];
        n = numel(platearea);
        G = [G; repmat(gnames(g), n,1)];
    end
    
    
    % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 
    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        % export anova
        vname = regexprep(t(1,:),'Prob>F','P_value');
        T = cell2table(t(2:end,:),'VariableNames',vname);
        cd(pSaveA); 
        writetable(T,[M{m},' Initial ANOVA.csv'],'Delimiter',',');
        
        % posthoc
        A.(M{m}).posthoc = multcompare_pairinterpretation(c,gnames);

        % export posthoc
        T = cell2table(a);
        cd(pSaveA); 
        writetable(T,[M{m},' Initial posthoc bonferroni.csv'],'Delimiter',',');
        
    end
end
MWTSet.Graph.HabRate_integral = A;


% Habituation rate (how fast does worms habituated to asymptote level; 
% 1 - asymptotic level is defined as no significant differences with the 
% responses to next 3 taps)
% 2- integral


% <200mM worms habituated slower (reversal responses duration and probability), higher doses increases rate (300mM+)



%% GRAPH: HABITUATION CURVES
pSaveA = [pSave,'/Graph Habituation curves'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
M = {'N_Total'
    'N_Reversed'
    'N_NoResponse'
    'RevFreq'
    'RevDur'
    'RevSpeed'};
Graph = MWTSet.Graph.HabCurve;
GroupName = regexprep(Graph.GroupNames,'_',' ');
timestamp = MWTSet.timestamp;
% get N2 
N2i = ~cellfun(@isempty,regexp(GroupName,'^N2'));

for m = 1:numel(M);% for each measure
    % get coordinates
    X = Graph.(M{m}).X;
    Y = Graph.(M{m}).Y;
    E = Graph.(M{m}).E;

    % create output table legend
    a = cellstr([repmat('t',size(Y,1),1) num2str([1:size(Y,1)]')]);
    a = regexprep(a,' ','');
    b = cellstr([repmat('SE',size(Y,1),1) num2str([1:size(Y,1)]')]);
    b = regexprep(b,' ','');
    vnames = [{'groupname'}; a;b];   
    % make table 
    d = [GroupName num2cell(Y') num2cell(E')];
    T = cell2table(d,'VariableNames',vnames);
    cd(pSaveA);
    writetable(T,[M{m},'.csv'],'Delimiter',',');
        
    % Create figure
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,'FontSize',18);
    box(axes1,'off');
    xlim(axes1,[0 size(Y,1)+1]);
    ylim(axes1,[min(min(Y-E))*.9 max(max(Y+E))*1.1]);
    hold(axes1,'all');
    errorbar1 = errorbar(X,Y,E,'LineWidth',1);
    color = [0 0 0; 1 0 0; [0.5 0.5 0.5]; [0.04 0.52 0.78]]; 
    n = min([size(color,1) size(Y,2)]);
    for x = 1:n
        set(errorbar1(x),'Color',color(x,1:3))
    end

    for x = 1:size(Y,2)
        set(errorbar1(x),'DisplayName',GroupName{x})
    end
    titlestr = timestamp;
    title(titlestr,'FontSize',12);
    xlabel('Tap','FontSize',18);
    ylabel(regexprep(M{m},'_',' '),'FontSize',18);
    legend1 = legend(axes1,'show');
    set(legend1,'EdgeColor',[1 1 1],...
        'Location','NorthEastOutside',...
        'YColor',[1 1 1],'XColor',[1 1 1],'FontSize',12);

    % text strings
    a = Graph.PlateN;
    b = '';
    for x = 1:numel(a)
        b = [b,num2str(a(x)),', '];
    end
    str1 = ['N(plate) = ',b(:,1:end-2)];
    textboxstr = [str1];
    annotation(figure1,'textbox',...
        [0.70 0.015 0.256 0.05],...
        'String',textboxstr,...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);
    % save fig
    savefigepsOnly(M{m},pSaveA);
    
end


%% GRAPH: HABITUATION CURVES BY STRAINS
pSaveA = [pSave,'/Graph Habituation curves by strain'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
M = {'N_Total','N_Reversed', 'N_NoResponse', 'RevFreq', 'RevDur','RevSpeed'};
Graph = MWTSet.Graph.HabCurve;
GroupName = regexprep(Graph.GroupNames,'_',' ');
timestamp = MWTSet.timestamp;
N2i = ~cellfun(@isempty,regexp(GroupName,'^N2')); % get N2 

if sum(N2i) ~= numel(GroupName)
    % get different strains
    a = regexpcellout(GroupName(~N2i),' ','split');
    strain = unique(a(:,1));

    for m = 1:numel(M);% for each measure
        % get coordinates
        X = Graph.(M{m}).X;
        Y = Graph.(M{m}).Y;
        E = Graph.(M{m}).E;

        % make figures by group and N2
        X1 = X; Y1 = Y;  E1 = E;
        for straini = 1:numel(strain)
            strainT = strain{straini};
            pSaveA2 = [pSaveA,'/',strainT];
            if isdir(pSaveA2) == 0; mkdir(pSaveA2); end
            iG = ~cellfun(@isempty,regexp(GroupName,['^',strainT]));
            gname = GroupName([find(N2i);find(iG)]);
            X = repmat(X1(:,1),1,sum(N2i)+sum(iG)); 
            Y = Y1(:,[find(N2i);find(iG)]); 
            E = E1(:,[find(N2i);find(iG)]);

            % create habituation curve
            figure1 = figure('Color',[1 1 1],'Visible','off');
            axes1 = axes('Parent',figure1,'FontSize',24);
            box(axes1,'off');
            xlim(axes1,[0 size(Y,1)+1]);
            ylim(axes1,[min(min(Y-E))*.9 max(max(Y+E))*1.1]);
            hold(axes1,'all');
            errorbar1 = errorbar(X,Y,E,'LineWidth',1);
            color = [0 0 0; 1 0 0; [0.5 0.5 0.5]; [0.04 0.52 0.78]]; 
            n = min([size(color,1) size(Y,2)]);
            for x = 1:n
                set(errorbar1(x),'Color',color(x,1:3))
            end

            for x = 1:size(Y,2)
                set(errorbar1(x),'DisplayName',gname{x})
            end
            titlestr = timestamp;
            title(titlestr,'FontSize',12);
            xlabel('Tap','FontSize',18);
            ylabel(regexprep(M{m},'_',' '),'FontSize',18);
            legend1 = legend(axes1,'show');
            set(legend1,'EdgeColor',[1 1 1],...
                'Location','NorthEastOutside',...
                'YColor',[1 1 1],'XColor',[1 1 1],'FontSize',12);

            % text strings
            a = Graph.PlateN(iG);
            b = '';
            for x = 1:numel(a)
                b = [b,num2str(a(x)),', '];
            end
            str1 = ['N(plate) = ',b(:,1:end-2)];
            textboxstr = [str1];
            annotation(figure1,'textbox',...
                [0.70 0.015 0.256 0.05],...
                'String',textboxstr,...
                'FitBoxToText','off',...
                'EdgeColor',[1 1 1]);
            % save fig
            savefigepsOnly150([M{m},' ',strain{straini}],pSaveA2);

        end
    end
end


%% GRAPH: HABITUATION CURVES BY STRAINS + STATS
pSaveA = [pSave,'/Graph HabCurve Combo by strain'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
M = {'RevFreq', 'RevDur','RevSpeed'};
Graph = MWTSet.Graph.HabCurve;
HabLevel = MWTSet.Graph.HabLevel;
Initial = MWTSet.Graph.Initial;
HabCurve = MWTSet.Graph.HabCurve;
HabRate = MWTSet.Graph.HabRate_integral;
GroupName = regexprep(Graph.GroupNames,'_',' ');
timestamp = MWTSet.timestamp;
N2i = ~cellfun(@isempty,regexp(GroupName,'^N2')); % get N2 
color = [0 0 0; 1 0 0; [0.04 0.52 0.78]; [0.47843137383461 0.062745101749897 0.894117653369904]]; 


if sum(N2i) ~= numel(GroupName)
    % get different strains
    a = regexpcellout(GroupName(~N2i),' ','split');
    strain = unique(a(:,1));

    for m = 1:numel(M);% for each measure
        % get coordinates
        X = Graph.(M{m}).X;
        Y = Graph.(M{m}).Y;
        E = Graph.(M{m}).E;

        % make figures by group and N2
        X1 = X; Y1 = Y;  E1 = E;
        for straini = 1:numel(strain)
            strainT = strain{straini};
            
            pSaveA2 = [pSaveA,'/',strainT];
            if isdir(pSaveA2) == 0; mkdir(pSaveA2); end
            
            %% create habituation curve
            iG = ~cellfun(@isempty,regexp(GroupName,['^',strainT]));
            gname = GroupName([find(N2i);find(iG)]);
            X = repmat(X1(:,1),1,sum(N2i)+sum(iG)); 
            Y = Y1(:,[find(N2i);find(iG)]); 
            E = E1(:,[find(N2i);find(iG)]);
            tapN = size(X,1);
            figure1 = figure('Color',[1 1 1],'Visible','off');
            axes1 = axes('Parent',figure1,'FontSize',24);
            box(axes1,'off');
            tapXlim = max(max(X))+1;
            xlim(axes1,[0 tapXlim+6]);
            
            rspMin = min(min(Y-E))*.9;
            rspMax = max(max(Y+E));
            ylim(axes1,[0 rspMax*1.1]);
            hold(axes1,'all');
            errorbar1 = errorbar(X,Y,E,'LineWidth',1,'Marker','o');
%             );
            n = min([size(color,1) size(Y,2)]);
            for x = 1:n
                set(errorbar1(x),'Color',color(x,1:3),...
                    'MarkerFaceColor',color(x,1:3),...
                'MarkerEdgeColor',color(x,1:3),...
                'LineStyle','none')
            end
            for x = 1:size(Y,2)
                set(errorbar1(x),'DisplayName',gname{x})
            end
            titlestr = timestamp;
            title(titlestr,'FontSize',12);
            xlabel('Tap','FontSize',18);
            ylabel(regexprep(M{m},'_',' '),'FontSize',18);
            legend1 = legend(axes1,'show');
            set(legend1,'EdgeColor',[1 1 1],...
                'Location','NorthEastOutside',...
                'YColor',[1 1 1],'XColor',[1 1 1],'FontSize',12);
            
            % text strings
            a = Graph.PlateN(iG);
            b = '';
            for x = 1:numel(a)
                b = [b,num2str(a(x)),', '];
            end
            str1 = ['N(plate) = ',b(:,1:end-2)];
            textboxstr = [str1];
            annotation(figure1,'textbox',...
                [0.70 0.015 0.256 0.05],...
                'String',textboxstr,...
                'FitBoxToText','off',...
                'EdgeColor',[1 1 1]);
            
            
            %% graph initial dot
            iG = ~cellfun(@isempty,regexp(MWTSet.Graph.Initial.(M{m}).GroupName,['^',strainT]));
            N2i = ~cellfun(@isempty,regexp(MWTSet.Graph.Initial.(M{m}).GroupName,'^N2')); % get N2
            % graph N2
            Y = []; E = []; X = [];
            Y(:,1) = Initial.(M{m}).Y(N2i);
            E(:,1) = Initial.(M{m}).E(N2i);
            X = repmat(tapXlim+1,size(Y));
            for x = 1:numel(Y)
                errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','o',...
                'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
                'MarkerFaceColor',color(x,1:3),...
                'MarkerEdgeColor',color(x,1:3));
            end
            % graph mutant
            Y = []; E = []; X = [];
            Y(:,1) = Initial.(M{m}).Y(iG);
            E(:,1) = Initial.(M{m}).E(iG);
            X = repmat(tapXlim+2,size(Y));
            for x = 1:numel(Y)
                errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','o',...
                'LineStyle','none','LineWidth',1,'Color',color(x+2,1:3),...
                'MarkerFaceColor',color(x+2,1:3),...
                'MarkerEdgeColor',color(x+2,1:3));
            end
            
            
            %% graph hab dot
            iG = ~cellfun(@isempty,regexp(Initial.(M{m}).GroupName,['^',strainT]));
            N2i = ~cellfun(@isempty,regexp(Initial.(M{m}).GroupName,'^N2')); % get N2
            % graph N2
            Y = []; E = []; X = [];
            Y(:,1) = HabLevel.(M{m}).Y(N2i);
            E(:,1) = HabLevel.(M{m}).E(N2i);
            X = repmat(tapXlim+1,size(Y));
            for x = 1:numel(Y)
                errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','v',...
                'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
                'MarkerFaceColor',color(x,1:3),...
                'MarkerEdgeColor',color(x,1:3));
            end
            % graph mutant
            Y = []; E = []; X = [];
            Y(:,1) = HabLevel.(M{m}).Y(iG);
            E(:,1) = HabLevel.(M{m}).E(iG);
            X = repmat(tapXlim+2,size(Y));
            for x = 1:numel(Y)
                errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','v',...
                'LineStyle','none','LineWidth',1,'Color',color(x+2,1:3),...
                'MarkerFaceColor',color(x+2,1:3),...
                'MarkerEdgeColor',color(x+2,1:3));
            end
            
 
            %% draw divider
            X = repmat(tapN+1,1,2);
            Y = [0 rspMax*1.1];
            line(X,Y,'Color',[0.3 0.3 0.3]);
            
            %% graph hab rate integral
            iG = ~cellfun(@isempty,regexp(HabRate.(M{m}).GroupName,['^',strainT]));
            N2i = ~cellfun(@isempty,regexp(HabRate.(M{m}).GroupName,'^N2')); % get N2
            idata = [find(N2i); find(iG)];
            % graph N2
            Y = []; E = []; 
            % Y will be expressed as average of area under the curve * 3 to
            % fit the graph display
            Y = HabRate.(M{m}).Y(N2i)./(tapN/1.5); 
            E = HabRate.(M{m}).E(N2i)./(tapN/1.5);
            X = repmat(tapXlim+4,size(Y));
            for x = 1:numel(Y)
                errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','x',...
                'LineStyle','none','LineWidth',1,'Color',color(x,1:3),...
                'MarkerFaceColor',color(x,1:3),...
                'MarkerEdgeColor',color(x,1:3));
            end
            % graph mutant
            Y = []; E = []; X = [];
            Y = HabRate.(M{m}).Y(iG)./(tapN/1.5); 
            E = HabRate.(M{m}).E(iG)./(tapN/1.5);
            X = repmat(tapXlim+5,size(Y));
            for x = 1:numel(Y)
                errorbar(X(x),Y(x),E(x),'MarkerSize',6,'Marker','x',...
                'LineStyle','none','LineWidth',1,'Color',color(x+2,1:3),...
                'MarkerFaceColor',color(x+2,1:3),...
                'MarkerEdgeColor',color(x+2,1:3));
            end
            
            %% draw divider
            X = repmat(tapN+4,1,2);
            Y = [0 rspMax*1.1];
            line(X,Y,'Color',[0.3 0.3 0.3]);
            % save fig
            savefigepsOnly150([M{m},' ',strain{straini}],pSaveA2);

        end
    end
end



%% EXCEL OUTPUT -GROUP
Data = MWTSet.Graph.HabCurve;
GroupNames = MWTSet.Graph.HabCurve.GroupNames;
MsrList = {'RevFreq','RevDur','RevSpeed'};
for Msri = 1:numel(MsrList)
    D = Data.(MsrList{Msri});
    d = [GroupNames num2cell(D.Y') num2cell(D.E')];

    a = cellstr([repmat('t',size(D.Y,1),1) num2str([1:size(D.Y,1)]')]);
    a = regexprep(a,' ','');
    b = cellstr([repmat('SE',size(D.Y,1),1) num2str([1:size(D.Y,1)]')]);
    b = regexprep(b,' ','');
    vnames = [{'groupname'}; a;b];
    
    T = cell2table(d,'VariableNames',vnames);
    cd(pSaveA);
    writetable(T,[MsrList{Msri},'.dat'],'Delimiter','\t');
end



%% GRAPH: INITIAL HAB LEVEL BAR
pSaveA = [pSave,'/','Graph Initial HabLevel bar'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end

gcategory = {'Initial','HabLevel'};


G = []; 
for m = 1:numel(M); % for each measure
    GroupName = {}; Y = []; E = []; X = [];
    for c = 1:numel(gcategory)
        GroupName(:,c) = MWTSet.Graph.(gcategory{c}).(M{m}).GroupName;
        X(:,c) = (1:numel(MWTSet.Graph.(gcategory{c}).(M{m}).Y))';
        Y(:,c) = MWTSet.Graph.(gcategory{c}).(M{m}).Y;
        E(:,c) = MWTSet.Graph.(gcategory{c}).(M{m}).E;
    end
    
    % create output table legend
    vnames = {'groupname', 'Initial','HabLevel','Initial_SE','HabLevel_SE'};   
    % make table 
    d = [GroupName(:,1) num2cell(Y) num2cell(E)];
    T = cell2table(d,'VariableNames',vnames);
    cd(pSaveA);
    writetable(T,[M{m},'Initial HabLevel.csv'],'Delimiter',',');    
  
    % GRAPHING
    titlename = MWTSet.timestamp;
    figname  = sprintf('%s Initial HabLevel bar',M{m});
%     figure1 = figure('Color',[1 1 1],'Visible','off');
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,...
        'XTickLabel',regexprep(GroupName(:,1),'_',' '),...
        'XTick',1:size(Y,1),...
        'FontSize',12);
    xlim(axes1,[0.5 size(Y,1)+0.5]);
    hold(axes1,'all');

    % create bar
    bar1 = bar(X,Y,'BarWidth',1,'Parent',axes1); 
    colorSet = [[0.043 0.52 0.78];[0.85 0.16 0]];
    for c = 1:size(colorSet,1)
        set(bar1(c),'FaceColor',colorSet(c,:),'DisplayName',gcategory{c});
    end

    % create errorbar
    errorbar1 = errorbar([X(:,1)-0.145 X(:,2)+0.145],Y,E);
    for c = 1:size(colorSet,1)
        set(errorbar1(c),...
            'LineStyle','none','Color',[0 0 0],'DisplayName',' ');
    end

    title(titlename);
    ylabel(M{m},'FontSize',16);
    legend1 = legend(axes1,'show');
    set(legend1,'EdgeColor',[1 1 1],...
        'Location','EastOutside','YColor',[1 1 1],'XColor',[1 1 1]);
    xticklabel_rotate; % rotate X labels
    savefigepsOnly(figname,pSaveA);    
end


%% GRAPH: INITAL VS HAB LEVEL DOT GRAPH
% if size(MWTSet.GraphSetting.GroupSeq,1) > 1

    pSaveA = [MWTSet.PATHS.pSaveA,'/','Graph Initial HabLevel dot'];
    if isdir(pSaveA) == 0; mkdir(pSaveA); end
    M = {'RevFreq', 'RevDur', 'RevSpeed'};
    for m = 1:numel(M)
        Y = []; E = []; X = [];
        gname = regexprep(MWTSet.Graph.Initial.(M{m}).GroupName,'_',' ');
        Y(:,1) = MWTSet.Graph.Initial.(M{m}).Y;
        E(:,1) = MWTSet.Graph.Initial.(M{m}).E;
        Y(:,2) = MWTSet.Graph.HabLevel.(M{m}).Y;
        E(:,2) = MWTSet.Graph.HabLevel.(M{m}).E;
        X(:,1:size(Y,2)) = repmat((1:size(Y,1))',1,size(Y,2));
    
        % make graph
        figure1 = figure('Color',[1 1 1],'Visible','off');
        axes1 = axes('Parent',figure1,...
            'XTickLabel',gname,...
            'XTick',1:size(Y,1),...
            'FontSize',14);
        hold(axes1,'all');
        errorbar1 = errorbar(X,Y,E,...
            'MarkerSize',6,'Marker','o',...
            'LineStyle','none',...
            'LineWidth',1);
        set(errorbar1(1),'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[0 0 0],...
            'DisplayName','Initial',...
            'Color',[0 0 0]);
        set(errorbar1(2),'MarkerFaceColor',[1 0 0],...
            'MarkerEdgeColor',[1 0 0],...
            'DisplayName','HabLevel',...
            'Color',[1 0 0]);
        ylabel(M{m},'FontSize',18);
        legend1 = legend(axes1,'show');
        set(legend1,'Location','EastOutside','EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);
%         savefigepsOnly150(['Initial vs hab level',M{m}],pSaveA);
        xticklabel_rotate; % rotate X labels
        savefigepsOnly([M{m},' Initial vs hab level'],pSaveA);    

    end
% end


%% GRAPH: HAB STRENGTH (INITIAL - HAB LEVEL)
pSaveA = [MWTSet.PATHS.pSaveA,'/','Graph Initial minus HabLevel'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
m = 1;
gnames = fieldnames(D);
for m = 1:numel(M)
    B = [];
    G = [];
    for g = 1:numel(gnames)
        % initial
        % hab level
        d = D.(gnames{g}).([M{m},'_Mean'])(end-2:end,:);
        hablevel = nanmean(d);
        initial = D.(gnames{g}).([M{m},'_Mean'])(1,:);
        d1 = initial - hablevel;
        % calculate last 3 taps average per plate 
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];
    end
    % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 

    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        % export anova
        vname = regexprep(t(1,:),'Prob>F','P_value');
        T = cell2table(t(2:end,:),'VariableNames',vname);
        cd(pSaveA); 
        writetable(T,[M{m},' Hab Strength ANOVA.csv'],'Delimiter',',');

        i = ((c(:,2) >0 & c(:,4) >0) + (c(:,2) <0 & c(:,4) <0)) >0;
        a = [gnames(c(:,1)), gnames(c(:,2))];
        a(i,3) = {'< 0.05'};
        a(~i,3) = {'not significant'};
        A.(M{m}).posthoc = a;
        % export posthoc
        T = cell2table(a);
        cd(pSaveA); 
        writetable(T,[M{m},' Hab Strength posthoc bonferroni.csv'],'Delimiter',',');
    end
end
MWTSet.Graph.HabStrength = A;

% graphing
Y = []; X = []; E = [];
for m = 1:numel(M)
    gname = regexprep(MWTSet.Graph.HabStrength.(M{m}).GroupName,'_',' ');
    Y = MWTSet.Graph.HabStrength.(M{m}).Y;
    E = MWTSet.Graph.HabStrength.(M{m}).E;
    X(:,1:size(Y,2)) = repmat((1:size(Y,1))',1,size(Y,2));
    
    
    % graph
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,...
        'XTickLabel',gname,...
        'XTick',1:size(Y,1),...
        'FontSize',14);
    hold(axes1,'all');
    errorbar1 = errorbar(X,Y,E,'MarkerSize',10,'Marker','o',...
        'LineStyle','none',...
        'LineWidth',1);
    set(errorbar1(1),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
        'DisplayName','Initial',...
        'Color',[0 0 0]);

    ylabel([M{m}, ' (initial - hab level)'],'FontSize',18);
%     legend1 = legend(axes1,'show');
%     set(legend1,'Location','EastOutside','EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);
%     savefigepsOnly150(['Habituation stength ',M{m}],pSaveA);
    xticklabel_rotate; % rotate X labels
    savefigepsOnly([M{m},' Hab strength'],pSaveA);  
end

%% GRAPH: HAB STRENGTH % (percent decrease Hab level - INITIAL)
pSaveA = [MWTSet.PATHS.pSaveA,'/','Graph Initial vs HabLevel'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
D = MWTSet.Data.ByGroupPerPlate;
m = 1;
gnames = fieldnames(D);
for m = 1:numel(M)
    B = [];
    G = [];
    for g = 1:numel(gnames)
        % hab level
        d = D.(gnames{g}).([M{m},'_Mean'])(end-2:end,:);
        hablevel = nanmean(d);
        % initial
        initial = D.(gnames{g}).([M{m},'_Mean'])(1,:);
        d1 = ((hablevel - initial)./initial)*100;
        % calculate last 3 taps average per plate 
        % get data
        B = [B;d1'];
        n = numel(d1);
        G = [G; repmat(gnames(g), n,1)];
        B(B == Inf) = nan;
    end
        % stats
    [m2, n2, se2,gnames2] = grpstats(B,G,{'mean','numel','sem','gname'});
    A.(M{m}).GroupName = gnames2;
    A.(M{m}).N = n2;
    A.(M{m}).Y = m2;
    A.(M{m}).E = se2; 

    if numel(unique(G)) > 1
        % anova
        [p,t,stats] = anova1(B,G,'off');
        [c,m1,h,gnames] = multcompare(stats,'ctype','bonferroni','display','off');
        A.(M{m}).ANOVA = t;
        A.(M{m}).posthoc = multcompare_pairinterpretation(c,gnames);

    end
end
MWTSet.Graph.HabStrength_percent = A;

% graphing
Y = []; X = []; E = [];

errorflag = 0;

for m = 1:numel(M)
    gnnow = MWTSet.Graph.HabStrength_percent.(M{m}).GroupName;
    if m > 1 && numel(gnnow) ~= size(Y,1)
        warning('some groups has missing %s data',M{m});
        errorflag = 1;
    else
        Y(:,m) = MWTSet.Graph.HabStrength_percent.(M{m}).Y;
        E(:,m) = MWTSet.Graph.HabStrength_percent.(M{m}).E;
        X(:,1:size(Y,2)) = repmat((1:size(Y,1))',1,size(Y,2));
    end
end    

gname = gnnow;

if errorflag == 0;
    figure1 = figure('Color',[1 1 1],'Visible','off');
    axes1 = axes('Parent',figure1,...
        'XTickLabel',regexprep(gname,'_',' '),...
        'XTick',1:size(Y,1),...
        'FontSize',12);
    hold(axes1,'all');
    errorbar1 = errorbar(X,Y,E,'MarkerSize',6,'Marker','o',...
        'LineStyle','none',...
        'LineWidth',1);
    colorset = [0 0 0; 1 0 0; 0 0 1];
    if size(colorset,1) > size(X,2)
        xN = size(X,2);
    else
        xN = size(colorset,1);
    end
    for x = 1:xN
        set(errorbar1(x),...
            'MarkerFaceColor',colorset(x,:),...
            'MarkerEdgeColor',colorset(x,:),...
            'DisplayName',M{x},...
            'Color',colorset(x,:));
    end

    ylabel('% diff (hab level - initial)','FontSize',16);
    legend1 = legend(axes1,'show');
    set(legend1,'Location','NorthOutside','EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);
%     savefigepsOnly150(['Habituation stength percent difference',M{m}],pSaveA);
    xticklabel_rotate; % rotate X labels
    savefigepsOnly('Hab stength percent difference',pSaveA);  
else
    warning('no [Habituation stength percent difference] graph created');
end



%% EXCEL OUTPUT
pSaveA = [MWTSet.PATHS.pSaveA,'/','Tables'];
if isdir(pSaveA) == 0; mkdir(pSaveA); end
gcategory = {'Initial','HabLevel'};
groupname = {};
rtype = {};
ctype = {};
N = []; Y = []; E = [];
for c = 1:numel(gcategory)
    for m = 1:numel(M)
        D = MWTSet.Graph.(gcategory{c}).(M{m});
        groupname = [groupname; D.GroupName];
        rtype = [rtype; repmat(M(m), numel(D.GroupName), 1)];
        ctype = [ctype; repmat(gcategory(c), numel(D.GroupName), 1)];
        N = [N; D.N];
        Y = [Y; D.Y];
        E = [E; D.E];       
    end
end
T = table;
T.group = groupname;
T.response_type = rtype;
T.msr = ctype;
T.N_plates = N;
T.Mean = Y;
T.SE = E;

cd(pSaveA);
writetable(T,'initial hab level summary.txt','Delimiter','\t');

% export anova 
for m = 1:numel(M)
    for c = 1:numel(gcategory)
        nn = fieldnames(MWTSet.Graph.(gcategory{c}).(M{m}));
        if strcmp(nn,'AVOVA') == 1
            D = MWTSet.Graph.(gcategory{c}).(M{m}).ANOVA;
            D(1,6) = {'p'};
            T = cell2table(D(2:end,:),'VariableNames',D(1,:));
            tname = sprintf('ANOVA %s %s.txt',gcategory{c},M{m});
            cd(pSaveA); 
            writetable(T, tname,'Delimiter','\t');
        end
    end   
end
% export posthoc
for m = 1:numel(M)
    for c = 1:numel(gcategory)
        nn = fieldnames(MWTSet.Graph.(gcategory{c}).(M{m}));
        if strcmp(nn,'posthoc') == 1
            D = MWTSet.Graph.(gcategory{c}).(M{m}).posthoc;
            T = cell2table(D,'VariableNames',{'group1','group2','p'});
            tname = sprintf('ANOVA %s %s posthoc bonferroni.txt',...
                gcategory{c},M{m});
            cd(pSaveA); 
            writetable(T, tname,'Delimiter','\t');
        end
    end   
end

%% SAVE MAT FILES
cd(pSave);
save('Dance_ShaneSpark3.mat','MWTSet');


%% Report done
fprintf('\n\n***DONE***\n\n');









