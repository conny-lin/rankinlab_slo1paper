function DataByPlate = initialResponseByPlate(pData)
%% calculate curve data by plate
%% Need to load data.mat, which includes MWTDB and DataMeta
load(pData);
pSave = fileparts(pData);
%% by plate
msru = {'curve','speedbm','speed'};
A = struct;
for msrui = 1:numel(msru)
    % get data
    D = DataMeta(:,{'mwtid',msru{msrui}});
    % get stats
    T = grpstats(D,'mwtid',{'numel','mean','sem'});
    % make table
    T1 = innerjoin(MWTDB(:,{'mwtid','groupname','expname','rx'}),T);
    % adjust table rx to dose
    i = ismember(T1.rx,'NA');
    T1.rx(i) = {'0mM'};
    % write table
    pF = fullfile(pSave,[msru{msrui},'_raw_byPlate.csv']);
    writetable(T1,pF);
    % put in struct
    A.(msru{msrui}) = T1;
   
    % get graph stats by plate
    x = T1(:,{'groupname',['mean_',msru{msrui}]});
    g = 'groupname';
    T = grpstats(x,g,{'numel','mean','sem'});
    % write table
    pF = fullfile(pSave,[msru{msrui},'_destats_byPlate.csv']);
    writetable(T,pF);
% %     [text,T,p,s,t,ST] = anova1_autoresults(T1,'groupname')
% %     [p,t,s] = anova1(T1.(['mean_',msru{msrui}]),T1.groupname,'off');
% %     text = anova_textresult(t);
%     x = T1.(['mean_',msru{msrui}]);
%     g = T1.groupname;
% %     [txt] = anova1_multcomp_auto(x,g)
%     pS = fullfile(pSave, [msru{msrui},'_stats_byPlate.txt']);
%     
%     anovan_std(x,g,{'strain','rx'},pS)
    return
end



DataByPlate = A;
save(pData,'DataByPlate', '-append');






















