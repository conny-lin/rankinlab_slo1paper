%% INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731
clc; clear; close all;
addpath('/Users/connylin/Dropbox/Code/Matlab/Library/General');
pSave = setup_std(mfilename('fullpath'),'RL','genSave',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731

%% GET DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731
pData = '/Users/connylin/Dropbox/Publication/Manuscript RL Alcohol hab model/Figures Tables Data/Fig6 response types/Fig7 60sISI response type/RspType_60sISI_20170703/RMANOVA.txt';
fileID = fopen(filename,'r');
D = textscan(fileID, '%s%[^\n\r]', 'Delimiter', '',  'ReturnOnError', false);
fclose(fileID);
D = [D{1:end-1}];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%20170731
