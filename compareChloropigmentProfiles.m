%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 15 15]);
set(0,'defaultAxesFontSize',10);

%% Load Data
cpigE = load("datafiles\chloro.mat","chloro300").chloro300;
cpigL = load('datafiles\chloro.mat','chloro300l').chloro300l;
pE = load('datafiles\chloro.mat','pres300').pres300(:,1);
pL = load('datafiles\chloro.mat','pres300l').pres300l(:,1);
time = load('datafiles\chloro.mat','days300').days300;

%% Calculate the mean profiles
cpigE_3m = mean(cpigE(:,1:30),2,"omitnan");
cpigE_300m = mean(cpigE(:,31:end),2,"omitnan");

%% Calculate the standard deviation

cpigE_3sd = std(cpigE(:,1:30),0,2,"omitnan");
cpigE_300sd = std(cpigE(:,31:end),0,2,"omitnan");

%% Plot

curve3a = cpigE_3m + cpigE_3sd;
curve3b = cpigE_3m - cpigE_3sd;
curve300a = cpigE_300m + cpigE_300sd;
curve300b = cpigE_300m - cpigE_300sd;

ax1 = figure;
plot(cpigE_3m,pE,'DisplayName','1988-1991 mean','Color','red','LineWidth',3);
hold on
plot(curve3a,pE,'DisplayName','1988-1991 STD','Color','red','LineStyle','--');
plot(curve3b,pE,'DisplayName','1988-1991 STD','Color','red','LineStyle','--','HandleVisibility','off');
plot(cpigE_300m,pE,'DisplayName','1991-2021 mean','Color','blue','LineWidth',3);
plot(curve300a,pE,'DisplayName','1988-1991 STD','Color','blue','LineStyle','--','HandleVisibility','off');
plot(curve300b,pE,'DisplayName','1988-1991 STD','Color','blue','LineStyle','--');
hold off
set(gca,'Ydir','reverse');
legend('Location','best');
xlabel('Chloropigment Concentration [{\mu}g L^{-1}]');
ylabel('Pressure [db]');
title('Comparison of Mean Concentration Profiles');
exportgraphics(ax1,'figures/comparisonOfProfiles.png')
clear ax1;