clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 28 12]);
set(0,'defaultAxesFontSize',10);

%% Load Test Data
chloro = load("datafiles\chloro.mat"','chloro256').chloro256;
chloroL = load("datafiles\chloro.mat",'chloro256_lang').chloro256_lang;
p = load('datafiles\chloro.mat','p_grid_256').p_grid_256(:,1);
p_lang = load('datafiles\chloro.mat','p_lang_grid256').p_lang_grid256(:,1);

time = load("datafiles\chloro.mat","time256").time256;

depth = linspace(0,-256,129);
%% median
medChloro = median(chloro,2,'omitnan');
medChloroL = median(chloroL,2,'omitnan');

%% Chloropigment SD
stdChloro = std(chloro,0,2,'omitnan');
stdChloroL = std(chloroL,0,2,'omitnan');

%% Chloropigment Skewness
skewChloro = skewness(chloro,0,2);
skewChloroL = skewness(chloroL,0,2);

%% Chloropigment Excess Kurtosis
kurtChloro = kurtosis(chloro,0,2);
kurtChloroL = kurtosis(chloroL,0,2);

%% Eulerian Moments

ax1 = figure;
subplot(1,4,1)
plot(medChloro,p);
set(gca, 'YDir','reverse');
xlim([0 1]); ylim([0 250]);
title('Chloropigments Median');
subplot(1,4,2)
plot(stdChloro,p);
set(gca, 'YDir','reverse');
xlim([0 0.2]); ylim([0 250]);
title('Chloropigments StD');
subplot(1,4,3)
plot(skewChloro,p);
set(gca, 'YDir','reverse');
xlim([0 4]); ylim([0 250]);
title('Chloropigments Skewness');
subplot(1,4,4)
plot(kurtChloro,p);
set(gca, 'YDir','reverse');
xlim([0 15]);
ylim([0 250]);
title('Chloropigments Excess Kurtosis');
exportgraphics(ax1,'figures/moments-Eulerian_1988_2021.png');

%% Lagrangian Moments

ax2 = figure;
subplot(1,4,1)
plot(medChloroL,p_lang);
set(gca, 'YDir','reverse');
ylim([-120 120]);
title('Chloropigments Median');
subplot(1,4,2)
plot(stdChloroL,p_lang);
set(gca, 'YDir','reverse');
ylim([-120 120]);
title('Chloropigments StD');
subplot(1,4,3)
plot(skewChloroL,p_lang);
set(gca, 'YDir','reverse');
ylim([-120 120]);
title('Chloropigments Skewness');
subplot(1,4,4)
plot(kurtChloroL,p_lang);
set(gca, 'YDir','reverse');
ylim([-120 120]);
title('Chloropigments Excess Kurtosis');
exportgraphics(ax2,'figures/moments-Lagrangian_1988_2021.png');
