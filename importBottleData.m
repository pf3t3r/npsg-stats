%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 15 15]);
set(0,'defaultAxesFontSize',10);

%% Load bottle data
chlBotData = importdata('data/hot-bottle.txt').data;
bottlePressure = chlBotData(:,4);
bottleChl = chlBotData(:,5);

%% Bin the chloropigment by pressure intervals

% bin such that central value of each bin is 5,10,15,etc.
binnedPressure = discretize(bottlePressure,2.5:5:202.5);

% Find no. of measurements at each depth
histSet = zeros(41,1);
x = 5:5:205;
for i = 1:41
    tmp = 0;
    for j = 1:6084
        if binnedPressure(j) == i
            tmp = tmp+1;
        end
    end
    histSet(i) = tmp;
    clear tmp;
end

figure;
barh(x,histSet);
set(gca,'YDir','reverse');