clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 28 12]);
set(0,'defaultAxesFontSize',10);

%% Load Test Data
chloro = load("datafiles\chloro.mat"','chloro256').chloro256;
chloroL = load("datafiles\chloro.mat",'chloro256l').chloro256l;
p = load('datafiles\chloro.mat','pgrid256').pgrid256(:,1);
p_lang = load('datafiles\chloro.mat','pgrid256l').pgrid256l(:,1);

time = load("datafiles\chloro.mat","time256").time256;

depth = linspace(0,-256,129);
n = length(depth);
%% Chloropigment Median
medChloro = median(chloro,2,'omitnan');
medChloroL = median(chloroL,2,'omitnan');

%% Chloropigment SD
stdChloro = std(chloro,0,2,'omitnan');
stdChloroL = std(chloroL,0,2,'omitnan');

% standard error of median = SD
y1_medE = medChloro + stdChloro;
y2_medE = medChloro - stdChloro;

% SE of SD
SE_stdC = stdChloro/(2*n - 2);
SE_stdCL = stdChloroL/(2*n - 2);

y1_stdE = stdChloro + SE_stdC;
y2_stdE = stdChloro - SE_stdC;

%% Chloropigment Skewness
skewChloro = skewness(chloro,0,2);
skewChloroL = skewness(chloroL,0,2);

% SE of Skewness: 
SES_C = sqrt((6*n*(n-1))/((n-2)*(n+1)*(n+3)));
SES_CL = SES_C;
y1_skeE = skewChloro + SES_C;
y2_skeE = skewChloro - SES_C;

%% Chloropigment Excess Kurtosis
kurtChloro = kurtosis(chloro,0,2);
kurtChloroL = kurtosis(chloroL,0,2);

% SE of Kurtosis
SEK_C = 4*(n^2-1)*SES_C / ((n-3)*(n+5));
y1_kurC = kurtChloro + SEK_C;
y2_kurC = kurtChloro - SEK_C;

%% Eulerian Moments

ax1 = figure;
subplot(1,4,1)
patch([y1_medE; flipud(y2_medE)], [p; flipud(p)], [0.8 0.8 0.8]);
hold on
plot(y1_medE,p,'b--');
plot(y2_medE,p,'b--');
plot(medChloro,p,'b-','LineWidth',2);
hold off
set(gca, 'YDir','reverse');
xlim([0 1]); ylim([0 250]);
title('Chloropigments Median');
subplot(1,4,2)
patch([y1_stdE; flipud(y2_stdE)], [p; flipud(p)], [0.8 0.8 0.8]);
hold on
plot(stdChloro,p,'b-','LineWidth',1);
plot(stdChloro + SE_stdC,p,'b--');
plot(stdChloro - SE_stdC,p,'b--');
hold off
set(gca, 'YDir','reverse');
xlim([0 0.2]); ylim([0 250]);
title('Chloropigments StD');
subplot(1,4,3)
patch([y1_skeE; flipud(y2_skeE)], [p; flipud(p)], [0.8 0.8 0.8]);
hold on
plot(skewChloro,p,'b-','LineWidth',2);
plot(skewChloro + SES_C,p,'b--');
plot(skewChloro - SES_C,p,'b--');
hold off
set(gca, 'YDir','reverse');
xlim([0 4]); ylim([0 250]);
title('Chloropigments Skewness');
subplot(1,4,4)
patch([y1_kurC-3; flipud(y2_kurC-3)], [p; flipud(p)], [0.8 0.8 0.8]);
hold on
plot(kurtChloro-3,p,'b-','LineWidth',2);
hold off
set(gca, 'YDir','reverse');
xlim([-1 15]);
ylim([0 250]);
title('Chloropigments Excess Kurtosis');
sgtitle('Chloropigment Moments (Eulerian)');
exportgraphics(ax1,'figures/moments-Eulerian_1988_2021.png');

%% Lagrangian Moments

ax2 = figure;
subplot(1,4,1)
plot(medChloroL,p_lang);
set(gca, 'YDir','reverse');
ylim([-120 120]); xlim([0 1]);
title('Chloropigments Median');
subplot(1,4,2)
plot(stdChloroL,p_lang);
set(gca, 'YDir','reverse');
ylim([-120 120]); xlim([0 0.2]);
title('Chloropigments StD');
subplot(1,4,3)
plot(skewChloroL,p_lang);
set(gca, 'YDir','reverse');
ylim([-120 120]); xlim([0 5]);
title('Chloropigments Skewness');
subplot(1,4,4)
plot(kurtChloroL-3,p_lang);
set(gca, 'YDir','reverse');
ylim([-120 120]); xlim([0 10]);
title('Chloropigments Excess Kurtosis');
sgtitle('Chloropigment Moments (Lagrangian)');
exportgraphics(ax2,'figures/moments-Lagrangian_1988_2021.png');

%% FILL example

x = [1 2 3 4 5]; % Both lines share same x value
y1 = x+1; % Equation for first line
y2 = 2*x; % Equation for second line

% plot the line edges
figure;
hold on 
plot(x, y1, 'LineWidth', 1);
plot(x, y2, 'LineWidth', 1);
% plot the shaded area
fill([x fliplr(x)], [y2 fliplr(y1)], 'r');

%% Another example
x = [1 2 3 3 2 2 1 1 2 2];
y = [10 10 20 30 40 50 50 40 30 20];

figure
fill(x,y,'r');