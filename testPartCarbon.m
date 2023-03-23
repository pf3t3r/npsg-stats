%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
set(0,'defaultAxesFontSize',10);

%% Load and Clean Bottle Data

p = importdata('data/particulateC.txt').data(:,4);
pC = importdata('data/particulateC.txt').data(:,5);
botid = importdata('data/particulateC.txt').data(:,1);

% Remove bottles taken at pressures below 2.5 db (within ~2.5 m of surface)
idRm = p > 2.5;
p = p(idRm);
pC = pC(idRm);
botid = botid(idRm);

% Remove bottles where particulate Carbon concentration <= 0
idZero = pC <= 0;
p = p(~idZero);
pC = pC(~idZero);
botid = botid(~idZero);

% Save cruise number (CRN) of each bottle - needed for Lagrangian analysis
tmp = num2str(botid);
bottleCRN = str2num(tmp(:,1:3));

clear tmp idRm idZero;

%% Bin Pressure

pBin = discretize(p,2.5:5:202.5);
% binRange = min(min(pBin)):10:max(max(pBin));
n = max(pBin);

%% Visualise
% figure
% plot(pC,p,'r.','DisplayName','Raw');
% hold on
% plot(pC,pBin,'b*','DisplayName','Binned');
% xlabel('Particulate Carbon [{\mu}mol kg^{-1}]');
% ylabel('Pressure [dbar]'); set(gca,'YDir','reverse');
% legend();

%% KS

ksEpC = nan(5,26);
obsPerLevel_E = nan(26,1);

for i = 1:n
    % find chl-a concentrations chla_i at binned pressure i
    pC_i = pC(pBin==i);
    % apply KS test to chla_i
    if length(pC_i) > 1
        [~,ksEpC(:,i),~] = statsplot2(pC_i,'noplot');
    end
    obsPerLevel_E(i) = length(pC_i);
    clear pC_i;
end
clear n;

%% KS FIG
depth5 = 5:5:200;

ax2 = figure;
subplot(1,2,1)
barh(obsPerLevel_E);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
ylim([0 40]);
set(gca,"YTick",1:1:40,"YTickLabel",depth5)
title('No. of Observations');

for i = 1:40
    if obsPerLevel_E(i) < 100
        ksEpC(:,i) = 0; % assign 0 instead of NaN b/c plot looks nicer :P
    end
end

subplot(1,2,2)
plot(ksEpC(1,:),depth5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEpC(2,:),depth5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpC(3,:),depth5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEpC(4,:),depth5,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Eulerian Bottle Particulate Carbon [5 dbar bins]');
exportgraphics(ax2,'figures/ks_bottle5dbar_particulateCarbon.png');
clear ax2;