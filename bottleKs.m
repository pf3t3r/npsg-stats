%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
set(0,'defaultAxesFontSize',10);

%% Load and Clean Bottle Data

bottlePressure = importdata('data/hotbot-88_21.txt').data(:,4);
bottleChl = importdata('data/hotbot-88_21.txt').data(:,5);
botid = importdata('data/hotbot-88_21.txt').data(:,1);

% Remove bottles taken at pressures below 2.5 db (within ~2.5 m of surface)
idRm = bottlePressure > 2.5;
bottlePressure = bottlePressure(idRm);
bottleChl = bottleChl(idRm);
botid = botid(idRm);

% Remove bottles where chl-a concentration = 0
idZero = bottleChl == 0;
bottlePressure = bottlePressure(~idZero);
bottleChl = bottleChl(~idZero);
botid = botid(~idZero);

% Save cruise number (CRN) of each bottle - needed for Lagrangian analysis
tmp = num2str(botid);
bottleCRN = str2num(tmp(:,1:3));
clear tmp;

% Remove bottles from cruises 330 on (b/c fluorescence analysis not done)
for i = 1:9929
    if bottleCRN(i) > 329
        id329 = i - 1;
        break;
    end
end

bottlePressure = bottlePressure(1:id329);
bottleChl = bottleChl(1:id329);
botid = botid(1:id329);

clear idRm idZero id329 i;
%% Bin Pressure

pBin5 = discretize(bottlePressure,2.5:5:202.5);
pBin10 = discretize(bottlePressure,2.5:10:202.5);
n = max(pBin5);
n10 = max(pBin10);

%% Visualise

% ax1 = figure;
% subplot(1,2,1)
% plot(bottleChl,bottlePressure,'r.','MarkerSize',4);
% xlabel('chl-a [ug/L]'); ylabel('Pressure [db]');
% set(gca,'YDir','reverse'); title('Raw Data');
% subplot(1,2,2)
% plot(bottleChl,pBin5,'b.','MarkerSize',4);
% xlabel('chl-a [ug/L]'); ylabel('Pressure [bin #]');
% set(gca,'YDir','reverse'); title('Binned Data');
% sgtitle('Fluorometric chl-a concentration (bottle data, 1988-2021)');
% exportgraphics(ax1,'figures/bottleChlaBinning.png');

%% Apply KS Test to chl-a across all pressures

ksEb5 = nan(5,40); ksEb10 = nan(5,20);
obsPerLevel_E5 = nan(40,1); obsPerLevel_E10 = nan(20,1);

for i = 1:n
    % find chl-a concentrations chla_i at binned pressure i
    chla_i = bottleChl(pBin5==i);
    % apply KS test to chla_i
    if length(chla_i) > 1
        [~,ksEb5(:,i),~] = statsplot2(chla_i,'noplot');
    end
    obsPerLevel_E5(i) = length(chla_i);
    clear chla_i;
end
clear n;

for i = 1:n10
    % find chl-a concentrations chla_i at binned pressure i
    chla_i10 = bottleChl(pBin10==i);
    % apply KS test to chla_i
    if length(chla_i10) > 1
        [~,ksEb10(:,i),~] = statsplot2(chla_i10,'noplot');
    end
    obsPerLevel_E10(i) = length(chla_i10);
    clear chla_i10;
end
clear n10;

depth5 = 5:5:200;
depth10 = 5:10:200;

ax2 = figure;
subplot(1,2,1)
barh(obsPerLevel_E5);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
ylim([0 40]);
set(gca,"YTick",1:1:40,"YTickLabel",depth5)
title('No. of Observations');

for i = 1:40
    if obsPerLevel_E5(i) < 100
        ksEb5(:,i) = 0; % assign 0 instead of NaN b/c plot looks nicer :P
    end
end
for i = 1:20
    if obsPerLevel_E10(i) < 100
        ksEb10(:,i) = 0; % assign 0 instead of NaN b/c plot looks nicer :P
    end
end


subplot(1,2,2)
plot(ksEb5(1,:),depth5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEb5(2,:),depth5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEb5(3,:),depth5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEb5(4,:),depth5,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Eulerian Bottle chl-a [5 dbar bins]');
exportgraphics(ax2,'figures/ks_bottle5dbar.png');
clear ax2;

%%

ax3 = figure;
subplot(1,2,1)
barh(obsPerLevel_E10);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
ylim([0.5 20]);
set(gca,"YTick",1:1:20,"YTickLabel",depth10)
title('No. of Observations');

subplot(1,2,2)
plot(ksEb10(1,:),depth10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEb10(2,:),depth10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEb10(3,:),depth10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEb10(4,:),depth10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Eulerian Bottle chl-a [10 dbar bins]');
exportgraphics(ax3,'figures/ks_bottle10dbar.png');
clear ax3;

%% Lagrangian Bottle

% apply adjustment according to cruise avg lagrangian adjustments
% need to import this adjustment 
% for this I need to save the isoL and ctdL structs as mat files in
% nightTimeKs.m

ctdL = load("datafiles\lagrangianData.mat").ctdL;

for i = 1:329
    tmp = ctdL(i).dcmID;
    tmp(tmp==1) = nan;
    dcmP(i) = 2*round(mean(tmp,'omitnan'));
end

% bottlePressure - dcmP = lagP
% max no. of bottles per crn ~30 (?)
tCRN = 1; tBotId = 1; % initialise with CRN1, with starting ID = 1
for i = 2:9792
    if bottleCRN(i) > bottleCRN(i-1)
        tCRN = [tCRN string(bottleCRN(i))]; % cruise no. 
        tBotId = [tBotId i]; % starting ID of bottles that cruise
    end
end


dcmP = dcmP(str2double(tCRN));

tChlaByCrn = nan(112,319);
tPByCrn = nan(112,319);
for i = 1:318
    tChlaByCrn(1:(tBotId(i+1)-tBotId(i)),i) = bottleChl(tBotId(i):(tBotId(i+1)-1));
    tPByCrn(1:(tBotId(i+1)-tBotId(i)),i) = bottlePressure(tBotId(i):(tBotId(i+1)-1));
end
tChlaByCrn(1:(9793-tBotId(319)),319) = bottleChl(tBotId(319):9792);
tPByCrn(1:(9793-tBotId(319)),319) = bottlePressure(tBotId(319):9792);

for i = 1:319
    pL(:,i) = tPByCrn(:,i) - dcmP(i);
end

pLBin10 = round(pL,-1);
% test3 = discretize(test1,-155:10:155);

% test3 = max(max(pLBin10));
% test4 = min(min(pLBin10));
pL_Range = min(min(pLBin10)):10:max(max(pLBin10));

clear bottleCRN ctdL tBotId tCRN;

%%

pL_levels = nan(627,31); % pressure
chla_levels = nan(627,31); % chla
obsPerLevel_L = [];

for i = 1:31
    tmp = []; tmp2 = [];
    % find where pLBin10 == -160
    for j = 1:112
        for k = 1:319
            if pLBin10(j,k) == pL_Range(i)
                tmp = [tmp pLBin10(j,k)];
                tmp2 = [tmp2 tChlaByCrn(j,k)];
            end
        end
    end
    pL_levels(1:length(tmp),i) = tmp';
    chla_levels(1:length(tmp),i) = tmp2';
    obsPerLevel_L = [obsPerLevel_L length(tmp)];
    clear tmp tmp2;
end
clear j k;

%% KS it all

ksLb10 = getKS(chla_levels,31);

for i = 1:31
    if obsPerLevel_L(i) < 100
        ksLb10(:,i) = 0; %again 0, not NaN, for nicer plot
    end
end
clear i;

%% figure this KS shit

ax4 = figure;

subplot(1,2,1)
histogram(pLBin10);
set(gca,'XDir','reverse');
set(gca,'view',[-90 90]);
set(gca,'XAxisLocation','top');
xlabel('Lagrangian Depth [dbar]');
title('No. of Observations')

subplot(1,2,2)
plot(ksLb10(1,:),pL_Range,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLb10(2,:),pL_Range,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLb10(3,:),pL_Range,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLb10(4,:),pL_Range,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle Chl-a');
exportgraphics(ax4,'figures/ksBotLag10db.png');
clear ax4;

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