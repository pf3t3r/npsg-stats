%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
set(0,'defaultAxesFontSize',10);

%% Load Bottle Data

bottlePressure = importdata('data/hotbot-88_21.txt').data(:,4);
bottleChl = importdata('data/hotbot-88_21.txt').data(:,5);

% To use only data from 1990 on, uncomment the following:
% bottlePressure = importdata('data\hotbot-90_21.txt').data(:,4);
% bottleChl = importdata('data\hotbot-90_21.txt').data(:,5);

%% Clean Data

% Remove bottles taken at pressures below 2.5 db (within ~2.5 m of surface)
idRm = bottlePressure > 2.5;
bottlePressure = bottlePressure(idRm);
bottleChl = bottleChl(idRm);

% Remove bottles where chl-a concentration = 0
idZero = bottleChl == 0;
bottlePressure = bottlePressure(~idZero);
bottleChl = bottleChl(~idZero);

clear idRm idZero;
%% Bin Pressure

binnedPressure = discretize(bottlePressure,2.5:5:202.5);
binnedPressure10 = discretize(bottlePressure,2.5:10:202.5);
n = max(binnedPressure);
n10 = max(binnedPressure10);

%% Visualise

ax1 = figure;
subplot(1,2,1)
plot(bottleChl,bottlePressure,'r.','MarkerSize',4);
xlabel('chl-a [ug/L]'); ylabel('Pressure [db]');
set(gca,'YDir','reverse'); title('Raw Data');
subplot(1,2,2)
plot(bottleChl,binnedPressure,'b.','MarkerSize',4);
xlabel('chl-a [ug/L]'); ylabel('Pressure [bin #]');
set(gca,'YDir','reverse'); title('Binned Data');
sgtitle('Fluorometric chl-a concentration (bottle data, 1988-2021)');
exportgraphics(ax1,'figures/bottleChlaBinning.png');

%% Apply KS Test to chl-a across all pressures

ksEb = nan(5,40); ksEb10 = nan(5,20);
measPerDepth = nan(40,1); measPerDepth10 = nan(40,1);

for i = 1:n
    % find chl-a concentrations chla_i at binned pressure i
    chla_i = bottleChl(binnedPressure==i);
    % apply KS test to chla_i
    if length(chla_i) > 1
        [~,ksEb(:,i),~] = statsplot2(chla_i);
    end
    measPerDepth(i) = length(chla_i);
    %clear chla_i;
end

for i = 1:n10
    % find chl-a concentrations chla_i at binned pressure i
    chla_i10 = bottleChl(binnedPressure10==i);
    % apply KS test to chla_i
    if length(chla_i10) > 1
        [~,ksEb10(:,i),~] = statsplot2(chla_i10);
    end
    measPerDepth10(i) = length(chla_i10);
    %clear chla_i;
end

x = 5:5:200;
x10 = 5:10:200;

ax2 = figure;
subplot(1,2,1)
barh(measPerDepth);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
ylim([0 40]);
set(gca,"YTick",1:1:40,"YTickLabel",x)
title('No. of Observations');

for i = 1:40
    if measPerDepth(i) < 100
        ksEb(:,i) = 0; % assign 0 instead of NaN b/c plot looks nicer :P
    end
end

subplot(1,2,2)
plot(ksEb(1,:),x,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEb(2,:),x,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEb(3,:),x,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEb(4,:),x,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian');

sgtitle('Kolmogorov Smirnov Test: Bottle chl-a [5 dbar bins]');
exportgraphics(ax2,'figures/ks_bottle5dbar.png');

%%

ax3 = figure;
subplot(1,2,1)
barh(measPerDepth10);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
ylim([0.5 20]);
set(gca,"YTick",1:1:20,"YTickLabel",x10)
title('No. of Observations');

subplot(1,2,2)
plot(ksEb10(1,:),x10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEb10(2,:),x10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEb10(3,:),x10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEb10(4,:),x10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian');

sgtitle('Kolmogorov Smirnov Test: Bottle chl-a [10 dbar bins]');
exportgraphics(ax3,'figures/ks_bottle10dbar.png');

%% Inspect for bimodality

% ax4 = figure;
% histfit(chla_i,[],'normal');
% hold on
% histfit(chla_i,[],'lognormal');
% histfit(chla_i,[],'weibull');
% histfit(chla_i,[],'gamma');
% legend();
% 
% exportgraphics(ax4,'figures/bottleLevelBimodality.png');