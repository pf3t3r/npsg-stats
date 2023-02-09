%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 15 15]);
set(0,'defaultAxesFontSize',10);

%% Load Data
chloro = load("datafiles\chloro.mat"','chloro256').chloro256;
chloroL = load("datafiles\chloro.mat"','chloro256_lang').chloro256_lang;
time = load("datafiles\chloro.mat","time256").time256;

%% Calculate the Kolmogorov-Smirnov Statistic for Chloropigment

depthMeasurements = 129;

mleE = zeros(5,2,depthMeasurements);
ksE = zeros(5,depthMeasurements);
nllE = zeros(5,depthMeasurements);

mleL = zeros(5,2,depthMeasurements);
ksL = zeros(5,depthMeasurements);
nllL = zeros(5,depthMeasurements);

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleE(:,:,i),ksE(:,i),nllE(:,i)] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleL(:,:,i),ksL(:,i),nllL(:,i)] = statsplot2(tmp,'noplot');
end

%% Plot the KS Statistic vs Depth for Five Distributions

ax1 = figure;
% Eulerian
subplot(1,2,1)
plot(ksE(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksE(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian KS-Test');

% Lagrangian
subplot(1,2,2)
plot(ksL(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksL(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian KS-Test');
exportgraphics(ax1,'figures/ks-EulAndLag_1988_2021.png');