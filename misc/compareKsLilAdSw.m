% Check effectiveness of four hypothesis data against ALOHA data.
% Compare a selection of hypothesis tests [1] against a selection of 
% parameters [2] from Station ALOHA.
% [1]: Kolmogorov-Smirnov (K-S), Lilliefors (Lil), Anderson-Darling (A-D),
% Shapiro-Wilks (S-W).
% [2]: Chl-a, Prochlorococcus, Synechococcus.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("packages\"); addpath("func\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

% Extract Maximum Mixed Layer Depth (per cruise) 'maxMld'
ctdData = importdata('datafiles\ctd_iso_ALL.mat').ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData i;
%% Chlorophyll a: 88-21

% 1. Load data
tmp = importdata('data/L1/hplcChla_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
ax = L1_compareHelper(tmp,maxMld);
sgtitle('[Chl a] 88-21: L1');
exportgraphics(ax,'figures/statTestComp/chla.png');
clear ax tmp;

%% Prochlorococcus: 05-21

% 1. Load data
tmp = importdata('data/L1/pro_05-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
ax2 = L1_compareHelper(tmp,maxMld,3);
sgtitle('Prochlorococcus 05-21: L1');
exportgraphics(ax2,'figures/statTestComp/pbact.png');
clear ax2 tmp;

%% Synechococcus: 05-21
% Two significant digits => may be unreliable!

% 1: Load data
tmp = importdata('data\L1\sbact_05-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
ax3 = L1_compareHelper(tmp,maxMld,3);
sgtitle('Synechococcus 05-21: L1');
exportgraphics(ax3,'figures/statTestComp/sbact.png');
clear ax3 tmp;