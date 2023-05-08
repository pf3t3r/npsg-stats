%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
set(0,'defaultAxesFontSize',10);

%% Basic Bin Depths

depth5 = 5:5:200;
depth10 = 5:10:200;

%% Load Bottle Data

% chla (regular method)
p_reg = importdata('data/hotbot-88_21.txt').data(:,4);
chl_reg = importdata('data/hotbot-88_21.txt').data(:,5);
id_reg = importdata('data/hotbot-88_21.txt').data(:,1);

% chla (HPLC method)
p_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,4);
chl_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,5);
id_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,1);

% ATP
p_atp = importdata('data/atp_88-21.txt').data(:,4);
atp = importdata('data/atp_88-21.txt').data(:,5);
id_atp = importdata('data/atp_88-21.txt').data(:,1);

%% Clean and Bin Data

[pBin5,pBin10,chl_Out,n5,n10] = cleanAndBin(p_reg,chl_reg,id_reg);                              % chla (regular method)
[pb5_hplc,pb10_hplc,chlOut_hplc,n5_hplc,n10_hplc] = cleanAndBin(p_hplc,chl_hplc,id_hplc);       % chla (HPLC method)
[pb5_atp,pb10_atp,atpOut,n5_atp,n10_atp] = cleanAndBin(p_atp,atp,id_atp);                       % ATP

%% Apply KS Test to chl-a across all pressures

[ks5, obs5, d5] = ksOfBinnedCon(chl_Out,pBin5,5);                    % 5 dbar / regular
[ks10, obs10, d10] = ksOfBinnedCon(chl_Out,pBin10,10);               % 10 dbar / regular
[ksHp5, obsHp5, dHp5] = ksOfBinnedCon(chlOut_hplc,pb5_hplc,5);       % 5 dbar / hplc
[ksHp10, obsHp10, dHp10] = ksOfBinnedCon(chlOut_hplc,pb10_hplc,10);  % 10 dbar / hplc

[ksAtp5, obsAtp5, dAtp5] = ksOfBinnedCon(atpOut,pb5_atp,5);          % 5 dbar / ATP
[ksAtp10, obsAtp10, dAtp10] = ksOfBinnedCon(atpOut,pb10_atp,10);     % 10 dbar / ATP

%% Eulerian 5 dbar: chl-a (regular method)

ax1 = figure;
subplot(1,2,1)
barh(obs5);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
set(gca,'YAxisLocation','right');
ylim([0 40]);
set(gca,"YTick",1:1:40,"YTickLabel",depth5)
title('No. of Observations');

subplot(1,2,2)
plot(ks5(1,:),d5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ks5(2,:),d5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks5(3,:),d5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks5(4,:),d5,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([0 200]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Regular Method: [chl a] (Eulerian Bottle, 5 dbar bin)');
exportgraphics(ax1,'figures/ks_bottleEulerian5db.png');
clear ax1;

%% Eulerian 5 dbar: chl-a (HPLC method)

ax2 = figure;
subplot(1,2,1)
barh(obsHp5);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
set(gca,'YAxisLocation','right');
ylim([0 40]);
set(gca,"YTick",1:1:40,"YTickLabel",depth5)
title('No. of Observations');

subplot(1,2,2)
plot(ksHp5(1,:),dHp5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksHp5(2,:),dHp5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksHp5(3,:),dHp5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksHp5(4,:),dHp5,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([0 200]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('HPLC Method: [chl a] (Eulerian Bottle, 5 dbar bin)');
exportgraphics(ax2,'figures/ks_HplcBottleEulerian5db.png');
clear ax2;

%% Eulerian 10 dbar: Regular Method

ax3 = figure;
subplot(1,2,1)
barh(obs10);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
set(gca,'YAxisLocation','right');
ylim([0.5 20]);
set(gca,"YTick",1:1:20,"YTickLabel",depth10)
title('No. of Observations');

subplot(1,2,2)
plot(ks10(1,:),d10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ks10(2,:),d10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks10(3,:),d10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks10(4,:),d10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Regular Method: [chl a] (Eulerian Bottle, 10 dbar bin)');
exportgraphics(ax3,'figures/ks_bottleEulerian10db.png');
clear ax3;


%% Eulerian 10 dbar: HPLC Method

ax4 = figure;
subplot(1,2,1)
barh(obsHp10);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
set(gca,'YAxisLocation','right');
ylim([0.5 20]);
set(gca,"YTick",1:1:20,"YTickLabel",depth10)
title('No. of Observations');

subplot(1,2,2)
plot(ksHp10(1,:),dHp10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksHp10(2,:),dHp10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksHp10(3,:),dHp10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksHp10(4,:),dHp10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('HPLC Method: [chl a] (Eulerian Bottle, 10 dbar bin)');
exportgraphics(ax4,'figures/ks_HplcBottleEulerian10db.png');
clear ax4;

%% Eulerian 5 dbar ATP

ax5 = figure;
subplot(1,2,1)
barh(obsAtp5);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
set(gca,'YAxisLocation','right');
ylim([0.5 40]);
set(gca,"YTick",1:1:40,"YTickLabel",depth5)
title('No. of Observations');

subplot(1,2,2)
plot(ksAtp5(1,:),dAtp5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksAtp5(2,:),dAtp5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksAtp5(3,:),dAtp5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksAtp5(4,:),dAtp5,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('ATP: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax5,'figures/ks_AtpEulerian5db.png');
clear ax5;

%% Eulerian 10 dbar ATP

ax6 = figure;
subplot(1,2,1)
barh(obsAtp10);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
set(gca,'YAxisLocation','right');
ylim([0.5 20]);
set(gca,"YTick",1:1:20,"YTickLabel",depth10)
title('No. of Observations');

subplot(1,2,2)
plot(ksAtp10(1,:),dAtp10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksAtp10(2,:),dAtp10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksAtp10(3,:),dAtp10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksAtp10(4,:),dAtp10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('ATP: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax6,'figures/ks_AtpEulerian10db.png');
clear ax6;

%% Hov -> works, just commented

% nb = 100;
% ax4 = figure;
% contourf(datenum(tTime),tPByCrn,tChlaByCrn,linspace(0,max(max(tChlaByCrn)),nb),'LineColor','auto');
% set(gca,'YDir','reverse');
% datetickzoom('x','yyyy mmm','keeplimits');
% colormap(flipud(cbrewer2('Spectral',nb)));
% c = colorbar;
% c.Label.String = 'Chlorophyll a [{\mu}g/L]';
% xlabel('Time');
% ylabel('Pressure [dbar]');
% title('Bottle Chlorophyll a: Eulerian (1988-2021)');
% 
% exportgraphics(ax4,'figures/hov_botChla.png');

%% Inspect for bimodality -> ehh

% ax4 = figure;
% histfit(chla_i,[],'normal');
% hold on
% histfit(chla_i,[],'lognormal');
% histfit(chla_i,[],'weibull');
% histfit(chla_i,[],'gamma');
% legend();
% 
% exportgraphics(ax4,'figures/bottleLevelBimodality.png');