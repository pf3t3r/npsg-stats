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
timeA = importdata('data/particulateC.txt').data(:,2); % mmddyy

% Remove bottles taken at pressures below 2.5 db (within ~2.5 m of surface)
idRm = p > 2.5;
p = p(idRm);
pC = pC(idRm);
botid = botid(idRm);
timeA = timeA(idRm);

% Remove bottles where particulate Carbon concentration <= 0
idZero = pC <= 0;
p = p(~idZero);
pC = pC(~idZero);
botid = botid(~idZero);
timeA = timeA(~idZero);

% Save cruise number (CRN) of each bottle - needed for Lagrangian analysis
tmp = num2str(botid);
bottleCRN = str2num(tmp(:,1:3));

clear tmp idRm idZero;

%% time transformation

timeA(timeA==-9) = 0;
tTest2 = num2str(timeA,'%06d');
tTest3 = tTest2(:,5:6); %year
tTest3a = str2num(tTest3);
tTest4 = tTest2(:,3:4); %dd
tTest5 = tTest2(:,1:2); %mm
tTest6 = [];
tTest6 = [tTest3 tTest5 tTest4];
tTest7 = str2num(tTest6);
tTest8 = datetime(tTest7,'ConvertFrom','yyyymmdd');

tTest9 = num2str(19*ones(965,1));
tTest10 = num2str(20*ones(2862-965,1));

tTest11 = [tTest9; tTest10];
tTest12 = [tTest11 tTest3];
tTest13 = [tTest12 tTest5 tTest4];
tTest14 = str2num(tTest13);
tTest14(tTest14==19000000) = NaN;
tTest14(tTest14==20000000) = NaN;
tTest15 = datetime(tTest14,'ConvertFrom','yyyymmdd');
tIme = tTest15;

%% Bin Pressure

pBin5 = discretize(p,2.5:5:202.5);
pBin10 = discretize(p,5:10:205);
n5 = max(pBin5); n10 = max(pBin10);

%% Visualise
% figure
% plot(pC,p,'r.','DisplayName','Raw');
% hold on
% plot(pC,pBin,'b*','DisplayName','Binned');
% xlabel('Particulate Carbon [{\mu}mol kg^{-1}]');
% ylabel('Pressure [dbar]'); set(gca,'YDir','reverse');
% legend();

%% Visualise Hovmoeller

% most frequent value in bottleCRN: used to define 2D Grid pressure
test1 = mode(bottleCRN);

% unique values in bottleCRN: used to define 2D Grid time
test2 = length(unique(bottleCRN)); 

% put time, pressure, and particulate carbon into 2D grids

%tCRN = 1;
tBotId = 1; % initialise with CRN1, with starting ID = 1
n = length(botid);
for i = 2:n
    if bottleCRN(i) > bottleCRN(i-1)
        %tCRN = [tCRN string(bottleCRN(i))]; % cruise no. 
        tBotId = [tBotId i]; % starting ID of bottles that cruise
    end
end

tpC = nan(test1,test2);
tP = nan(test1,test2);
tTime = NaT(test1,test2);
for i = 1:test2-1
    tpC(1:(tBotId(i+1)-tBotId(i)),i) = pC(tBotId(i):(tBotId(i+1)-1));
    tP(1:(tBotId(i+1)-tBotId(i)),i) = p(tBotId(i):(tBotId(i+1)-1));
    tTime(1:(tBotId(i+1)-tBotId(i)),i) = tIme(tBotId(i):(tBotId(i+1)-1));
end
tpC(1:(n+1-tBotId(309)),319) = pC(tBotId(309):n);
tP(1:(n+1-tBotId(309)),319) = p(tBotId(309):n);
tTime(1:(n+1-tBotId(309)),319) = tIme(tBotId(309):n);

clear n;

nb = 100;
ax1 = figure;
contourf(datenum(tTime),tP,tpC,linspace(0,max(max(tpC)),nb),'LineColor','auto');
set(gca,'YDir','reverse');
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Particulate Carbon [{\mu}mol/kg]';
xlabel('Time');
ylabel('Pressure [dbar]');
title('Bottle Particulate Carbon: Eulerian (1989-2020)');

exportgraphics(ax1,'figures/hov_partC.png');
%% KS

% 5 DBAR BINS
ksEpC5 = nan(5,40);
obsPerLevel_E5 = nan(40,1);

for i = 1:n5
    % find chl-a concentrations chla_i at binned pressure i
    pC_i = pC(pBin5==i);
    % apply KS test to chla_i
    if length(pC_i) > 1
        [~,ksEpC5(:,i),~] = statsplot2(pC_i,'noplot');
    end
    obsPerLevel_E5(i) = length(pC_i);
    clear pC_i;
end
clear n5;

ksEpC10 = nan(5,21);
obsPerLevel_E10 = nan(21,1);

for i = 1:n10
    % find chl-a concentrations chla_i at binned pressure i
    pC_i = pC(pBin10==i);
    % apply KS test to chla_i
    if length(pC_i) > 1
        [~,ksEpC10(:,i),~] = statsplot2(pC_i,'noplot');
    end
    obsPerLevel_E10(i) = length(pC_i);
    clear pC_i;
end
clear n10;
%% FIG: Eulerian KS 5 dbar
depth5 = 5:5:200; d5line = depth5;

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
        ksEpC5(:,i) = nan; % assign nan
        d5line(i) = nan;
    end
end

subplot(1,2,2)
plot(ksEpC5(1,:),depth5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
yline(d5line,'Color',[0.6 0.6 0.6],'HandleVisibility','off');
plot(ksEpC5(2,:),depth5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpC5(3,:),depth5,'xr-','DisplayName','Weibull','MarkerSize',6);
plot(ksEpC5(4,:),depth5,'r.--','DisplayName','Gamma','MarkerSize',6);
plot(max(ksEpC5),depth5,'gO','DisplayName','Most Likely','MarkerSize',8);
hold off
grid off;
ax = gca; ax.XGrid = 'on';
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Eulerian Bottle Particulate Carbon [5 dbar bins]');
exportgraphics(ax2,'figures/ks_bottle5dbar_particulateCarbon.png');
clear ax2;

%% FIG: Eulerian KS 10 dbar

depth10 = 5:10:205; d10line = depth10;

ax3 = figure;
subplot(1,2,1)
barh(obsPerLevel_E10);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
ylim([0 20]);
set(gca,"YTick",1:1:20,"YTickLabel",depth10)
title('No. of Observations');

for i = 1:21
    if obsPerLevel_E10(i) < 100
        ksEpC10(:,i) = nan;
        d10line(i) = nan;
    end
end

subplot(1,2,2)
plot(ksEpC10(1,:),depth10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
yline(d5line,'Color',[0.6 0.6 0.6],'HandleVisibility','off');
plot(ksEpC10(2,:),depth10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpC10(3,:),depth10,'xr-','DisplayName','Weibull','MarkerSize',6);
plot(ksEpC10(4,:),depth10,'r.--','DisplayName','Gamma','MarkerSize',6);
plot(max(ksEpC10),depth10,'gO','DisplayName','Most Likely','MarkerSize',8);
hold off
grid off;
ax = gca; ax.XGrid = 'on';
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Eulerian Bottle Particulate Carbon [10 dbar bin]');
exportgraphics(ax3,'figures/ks_bottle10dbar_particulateCarbon.png');
clear ax3;