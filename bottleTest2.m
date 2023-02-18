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
n = max(binnedPressure);

%% Visualise

figure;
subplot(1,2,1)
plot(bottleChl,bottlePressure,'r.','MarkerSize',4);
xlabel('chl-a [ug/L]'); ylabel('Pressure [db]');
set(gca,'YDir','reverse'); title('Raw Data');
subplot(1,2,2)
plot(bottleChl,binnedPressure,'b.','MarkerSize',4);
xlabel('chl-a [ug/L]'); ylabel('Pressure [bin #]');
set(gca,'YDir','reverse'); title('Binned Data');
sgtitle('Fluorometric chl-a concentration (bottle data, 1988-2021)');

%% Apply KS Test to chl-a across all pressures

ksEb = nan(5,40);
measPerDepth = [];
for i = 1:n
    % find chl-a concentrations chla_i at binned pressure i
    chla_i = bottleChl(binnedPressure==i);
    % apply KS test to chla_i
    if length(chla_i) > 1
        [~,ksEb(:,i),~] = statsplot2(chla_i,'noplot');
    end
    measPerDepth = [measPerDepth length(chla_i)];
    clear chla_i;
end

x = linspace(5,200,n);

figure;
subplot(1,2,1)
barh(measPerDepth);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
ylim([0 40]);
set(gca,"YTick",1:1:40,"YTickLabel",x)
title('No. of Observations');

subplot(1,2,2)
plot(ksEb(1,:),x,'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksEb(2,:),x,'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksEb(3,:),x,'Color','red','DisplayName','Weibull');
plot(ksEb(4,:),x,'Color','red','LineStyle','--','DisplayName','Gamma');
hold off
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian KS-Test');