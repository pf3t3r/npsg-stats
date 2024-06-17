% Draw mathematical moments for chloropigment time series.
% Moments: mean, standard deviation, skewness, and (excess) kurtosis.

clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on'); set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 28 12]);
set(0,'defaultAxesFontSize',10);

%% Load Test Data
chlE = load("datafiles\chloro.mat"','chloro256').chloro256;
% chlL = load("datafiles\chloro.mat",'chloro256l').chloro256l;
pE = load('datafiles\chloro.mat','pgrid256').pgrid256(:,1);
% pL = load('datafiles\chloro.mat','pgrid256l').pgrid256l(:,1);

n = length(pE);

%% Chloropigment Median: E+L
medChloro = median(chlE,2,'omitnan');
% medChloroL = median(chlL,2,'omitnan');

%% Chloropigment SD

% EULERIAN
% Standard Deviation
stdChloro = std(chlE,0,2,'omitnan');
% Standard Error of Median (= Standard Deviation)
y1_medE = medChloro + stdChloro;
y2_medE = medChloro - stdChloro;
% Standard Error of Standard Deviation
SE_stdE = stdChloro/(2*n - 2);
y1_stdE = stdChloro + SE_stdE;
y2_stdE = stdChloro - SE_stdE;

% % LAGRANGIAN
% % Standard Deviation
% stdChloroL = std(chlL,0,2,'omitnan');
% % Standard Error of Median (= Standard Deviation)
% y1_medL = medChloroL + stdChloroL;
% y2_medL = medChloroL - stdChloroL;
% % Standard Error of Standard Deviation
% SE_stdL = stdChloroL/(2*n - 2);
% y1_stdL = stdChloroL + SE_stdL;
% y2_stdL = stdChloroL - SE_stdL;

%% Chloropigment Skewness

% EULERIAN
% Skewness
skewChloro = skewness(chlE,0,2);
% Standard Error of Skewness: 
SES_E = sqrt((6*n*(n-1))/((n-2)*(n+1)*(n+3)));
y1_skeE = skewChloro + SES_E;
y2_skeE = skewChloro - SES_E;

% % LAGRANGIAN
% % Skewness
% skewChloroL = skewness(chlL,0,2);
% % Standard Error of Skewness
% SES_L = SES_E;
% y1_skeL = skewChloroL + SES_L;
% y2_skeL = skewChloroL - SES_L;

%% Chloropigment Excess Kurtosis

% EULERIAN
% Excess Kurtosis
kurtChloro = kurtosis(chlE,0,2);
% Standard Error on Excess Kurtosis
SEK_E = 4*(n^2-1)*SES_E / ((n-3)*(n+5));
y1_kurE = kurtChloro + SEK_E;
y2_kurE = kurtChloro - SEK_E;

% % LAGRANGIAN
% % Excess Kurtosis
% kurtChloroL = kurtosis(chlL,0,2);
% % Standard Error on Excess Kurtosis
% SEK_L = SEK_E;
% y1_kurL = kurtChloroL + SEK_L;
% y2_kurL = kurtChloroL - SEK_L;


%% Average Profile

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 10 12]);

ax0 = figure;
patch([y1_medE; flipud(y2_medE)], [pE; flipud(pE)], [0.8 0.8 0.8]);
hold on
plot(medChloro,pE,'LineWidth',1);
hold off
set(gca, 'YDir','reverse');
% xlim([0 1]);
ylim([0 250]);
ylabel('Pressure [dbar]',FontSize=15); xlabel('chloropigment [$\mu$g/L]',Interpreter='latex',FontSize=15);
% title('Median');

exportgraphics(ax0,'figures/moments-Eulerian_1988_2021.png');

%% Eulerian Moments
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 28 12]);

ax1 = figure;
sgtitle('Chloropigment Moments (Eulerian)');

subplot(1,4,1)
patch([y1_medE; flipud(y2_medE)], [pE; flipud(pE)], [0.8 0.8 0.8]);
hold on
plot(medChloro,pE,'LineWidth',1);
hold off
set(gca, 'YDir','reverse');
xlim([0 1]); ylim([0 250]);
title('Median');

subplot(1,4,2)
patch([y1_stdE; flipud(y2_stdE)], [pE; flipud(pE)], [0.8 0.8 0.8]);
hold on
plot(stdChloro,pE,'LineWidth',1);
hold off
set(gca, 'YDir','reverse');
xlim([0 0.2]); ylim([0 250]);
title('StD');

subplot(1,4,3)
patch([y1_skeE; flipud(y2_skeE)], [pE; flipud(pE)], [0.8 0.8 0.8]);
hold on
plot(skewChloro,pE,'LineWidth',1);
hold off
set(gca, 'YDir','reverse');
xlim([0 4]); ylim([0 250]);
title('Skewness');

subplot(1,4,4)
patch([y1_kurE-3; flipud(y2_kurE-3)], [pE; flipud(pE)], [0.8 0.8 0.8]);
hold on
plot(kurtChloro-3,pE,'LineWidth',1);
hold off
set(gca, 'YDir','reverse');
xlim([-1 15]);
ylim([0 250]);
title('Excess Kurtosis');

exportgraphics(ax1,'figures/moments-Eulerian_1988_2021.png');

% %% Lagrangian Moments
% 
% ax2 = figure;
% sgtitle('Chloropigment Moments (Lagrangian)');
% 
% subplot(1,4,1)
% patch([y1_medL; flipud(y2_medL)], [pL; flipud(pL)], [0.8 0.8 0.8]);
% hold on
% plot(medChloroL,pL,'LineWidth',1);
% hold off
% set(gca, 'YDir','reverse');
% ylim([-120 120]); xlim([0 1]);
% title('Median');
% 
% subplot(1,4,2)
% patch([y1_stdL; flipud(y2_stdL)], [pL; flipud(pL)], [0.8 0.8 0.8]);
% hold on
% plot(stdChloroL,pL,'LineWidth',1);
% hold off
% set(gca, 'YDir','reverse');
% ylim([-120 120]); xlim([0 0.2]);
% title('StD');
% 
% subplot(1,4,3)
% patch([y1_skeL; flipud(y2_skeL)], [pL; flipud(pL)], [0.8 0.8 0.8]);
% hold on
% plot(skewChloroL,pL,'LineWidth',1);
% hold off
% set(gca, 'YDir','reverse');
% ylim([-120 120]); xlim([0 5]);
% title('Skewness');
% 
% subplot(1,4,4)
% patch([y1_kurL-3; flipud(y2_kurL-3)], [pL; flipud(pL)], [0.8 0.8 0.8]);
% hold on
% plot(kurtChloroL-3,pL,'LineWidth',1);
% hold off
% set(gca, 'YDir','reverse');
% ylim([-120 120]); xlim([0 10]);
% title('Excess Kurtosis');
% 
% exportgraphics(ax2,'figures/moments-Lagrangian_1988_2021.png');