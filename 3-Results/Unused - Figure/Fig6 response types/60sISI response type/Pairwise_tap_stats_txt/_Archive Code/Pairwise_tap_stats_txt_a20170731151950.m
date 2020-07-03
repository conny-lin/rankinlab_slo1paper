%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731

%% GET DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731
pData = '/Users/connylin/Dropbox/Publication/Manuscript RL Alcohol hab model/Figures Tables Data/Fig6 response types/Fig7 60sISI response type/RspType_60sISI_20170703/RMANOVA.txt';
fileID = fopen(pData,'r');
D = textscan(fileID, '%s%[^\n\r]', 'Delimiter', '',  'ReturnOnError', false);
fclose(fileID);
D = [D{1:end-1}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731

%% GET TXT FROM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731
% get index
i = regexpcellout(D,'^[*]{3}\s{1}Posthoc[(]Tukey[)]group[*]response_type\s{1}by\s{1}tap\s{1}[*]{3}');
i = find(i);
i = i(end);
% get only data after index
D = D(i:end);
%% split info by *
A = regexpcellout(D,'(*)|(,\sp\s)','split');
% get tap
taps = A(:,1);
taps = regexprep(taps,'tap','');
% get group 1 group
g1 = A(:,2);
% get group 1 type
t1 = A(:,3);
% get group 2 group
g2 = A(:,4);
% get group 2 type
t2 = A(:,5);
% get p value
p = A(:,6);

%% group by p value
% not significant
pns = regexpcellout(p,'n[.]s[.]');
% less than 0.001
p001 = regexpcellout(p,'<\s[.]001');
% = 0.001 or others
peq = regexpcellout(p,'=\s[.]');

%% get groups by signiifcant and not sig
i = pns;
PNS = [taps(i) g1(i) t1(i) g2(i) t2(i) p(i)];
i = p001 | peq;
PS = [taps(i) g1(i) t1(i) g2(i) t2(i) p(i)];

% create list of taps
tlist = 1:30;
tlist = num2cellstr(tlist');

%% get rev vs acc, N2 vs 400mM
g1t = 'N2';
g2t = 'N2';
t1t = 'rev';
t2t = 'acc';

B = PS;
i = ismember(B(:,2),g1t) & ismember(B(:,4),g2t) & ismember(B(:,3),t1t) & ismember(B(:,5),t2t);
j = ismember(B(:,2),g1t) & ismember(B(:,4),g2t) & ismember(B(:,3),t2t) & ismember(B(:,5),t1t);
R = B(i|j,:);
% check p values unique
a = setdiff(R(:,1),tlist);
if numel(unique(R(:,6))) == 1
    ptxt = unique(R(:,6));
    if isempty(a)
        taptxt = {'all'};
        rtext = strjoin([taptxt ptxt],', p');
    end
else % go by each unique value
    
end
% find taps
rtime = cellfun(@str2num,R(:,1));
a = R(:,6);
a(ismember(R(:,6),'< .001')) = {'= .0001'};
pv = cellfun(@str2num,regexprep(a,'((<)|(=))\s',''));
[pstr,pvt,tv] = pvaluestring(rtime,pv,0.05,0.001,'space',true,'zerob4',false);
pstr0 = pstr;



%% get rev vs acc, N2 vs 400mM
g1t = 'N2_400mM';
g2t = 'N2_400mM';
t1t = 'rev';
t2t = 'acc';

B = PS;
i = ismember(B(:,2),g1t) & ismember(B(:,4),g2t) & ismember(B(:,3),t1t) & ismember(B(:,5),t2t);
j = ismember(B(:,2),g1t) & ismember(B(:,4),g2t) & ismember(B(:,3),t2t) & ismember(B(:,5),t1t);
R = B(i|j,:);
% check p values unique
a = setdiff(R(:,1),tlist);
if numel(unique(R(:,6))) == 1
    ptxt = unique(R(:,6));
    if isempty(a)
        taptxt = {'all'};
        rtext = strjoin([taptxt ptxt],', p');
    end
else % go by each unique value
    
end
% find taps
rtime = cellfun(@str2num,R(:,1));
a = R(:,6);
a(ismember(R(:,6),'< .001')) = {'= .0001'};
pv = cellfun(@str2num,regexprep(a,'((<)|(=))\s',''));
[pstr,pvt,tv] = pvaluestring(rtime,pv,0.05,0.001,'space',true,'zerob4',false);
pstr4 = pstr;

%% summary
fid = fopen(fullfile(pSave,'result.txt'),'w');
fprintf('0mM\n%s 400mM\n%s\n ',pstr0,pstr4);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731























