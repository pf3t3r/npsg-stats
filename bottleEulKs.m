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

% chla monovinyl
id_cmo = importdata('data/chlaMonovinyl_88-21.txt').data(:,1);
p_cmo = importdata('data/chlaMonovinyl_88-21.txt').data(:,4);
cmo = importdata('data/chlaMonovinyl_88-21.txt').data(:,5);

% HPLC Chlorophyll a Divinyl
id_cdi = importdata('data/chlaDivinyl_88-21.txt').data(:,1);
p_cdi = importdata('data/chlaDivinyl_88-21.txt').data(:,4);
cdi = importdata('data/chlaDivinyl_88-21.txt').data(:,5);

% ATP
p_atp = importdata('data/atp_88-21.txt').data(:,4);
atp = importdata('data/atp_88-21.txt').data(:,5);
id_atp = importdata('data/atp_88-21.txt').data(:,1);

% Particulate Carbon
p_parc = importdata('data/parC_89-20.txt').data(:,4);
parc = importdata('data/parC_89-20.txt').data(:,5);
id_parc = importdata('data/parC_89-20.txt').data(:,1);

% Particulate Nitrogen
p_nit = importdata('data/parN_89-21.txt').data(:,4);
nit = importdata('data/parN_89-21.txt').data(:,5);
id_nit = importdata('data/parN_89-21.txt').data(:,1);

% Particulate Phosphorus
p_pho = importdata('data/parP_88-21.txt').data(:,4);
pho = importdata('data/parP_88-21.txt').data(:,5);
id_pho = importdata('data/parP_88-21.txt').data(:,1);

% Particulate Phosphorus - after 2011
p_pho11 = importdata('data/parP_11-21.txt').data(:,4);
pho11 = importdata('data/parP_11-21.txt').data(:,5);
id_pho11 = importdata('data/parP_11-21.txt').data(:,1);

% Heterotrophic Bacteria
id_het = importdata('data/hetBac_05-21.txt').data(:,1);
p_het = importdata('data/hetBac_05-21.txt').data(:,4);
het = importdata('data/hetBac_05-21.txt').data(:,5);

% Prochlorococcus
id_pro = importdata('data/prochl_05-21.txt').data(:,1);
p_pro = importdata('data/prochl_05-21.txt').data(:,4);
pro = importdata('data/prochl_05-21.txt').data(:,5);

% Synechococcus
id_syn = importdata('data/synech_05-21.txt').data(:,1);
p_syn = importdata('data/synech_05-21.txt').data(:,4);
syn = importdata('data/synech_05-21.txt').data(:,5);

% Pico-Eukaryotes
id_pic = importdata('data/picoeu_05-21.txt').data(:,1);
p_pic = importdata('data/picoeu_05-21.txt').data(:,4);
pic = importdata('data/picoeu_05-21.txt').data(:,5);

%% Clean and Bin Data

[pBin5,pBin10,chl_Out,n5,n10] = cleanAndBin(p_reg,chl_reg,id_reg);                              % Chlorophyll a (regular method)
[pb5_hplc,pb10_hplc,chlOut_hplc,n5_hplc,n10_hplc] = cleanAndBin(p_hplc,chl_hplc,id_hplc);       % Chlorophyll a (HPLC method)
[pb5_atp,pb10_atp,atpOut,n5_atp,n10_atp] = cleanAndBin(p_atp,atp,id_atp);                       % ATP
[pb5_cmo,pb10_cmo,cmoOut,n5_cmo,n10_cmo] = cleanAndBin(p_cmo,cmo,id_cmo);                       % HPLC Monovinyl Chlorophyll a
[pb5_cdi,pb10_cdi,cdiOut,n5_cdi,n10_cdi] = cleanAndBin(p_cdi,cdi,id_cdi);                       % HPLC Divinyl Chlorophyll a
[pb5_parc,pb10_parc,parcOut,n5_parc,n10_parc] = cleanAndBin(p_parc,parc,id_parc);               % Particulate Carbon
[pb5_nit,pb10_nit,nitOut,n5_nit,n10_nit] = cleanAndBin(p_nit,nit,id_nit);                       % Particulate Nitrogen
[pb5_pho,pb10_pho,phoOut,n5_pho,n10_pho] = cleanAndBin(p_pho,pho,id_pho);                       % Particulate Phosphorus
[pb5_pho11,pb10_pho11,phoOut11,n5_pho11,n10_pho11] = cleanAndBin(p_pho11,pho11,id_pho11);                       % Particulate Phosphorus (11-21)
[pb5_het,pb10_het,hetOut,n5_het,n10_het] = cleanAndBin(p_het,het,id_het);                       % Heterotrophic Bacteria
[pb5_pro,pb10_pro,proOut,n5_pro,n10_pro] = cleanAndBin(p_pro,pro,id_pro);                       % Prochlorococcus
[pb5_syn,pb10_syn,synOut,n5_syn,n10_syn] = cleanAndBin(p_syn,syn,id_syn);                       % Synechococcus
[pb5_pic,pb10_pic,picOut,n5_pic,n10_pic] = cleanAndBin(p_pic,pic,id_pic);                       % Pico-Eukaryotes

%% Apply KS Test to chl-a across all pressures

[ks5, obs5, d5] = ksOfBinnedCon(chl_Out,pBin5,5);                    % 5 dbar / Chlorophyll a (Regular Method)
[ks10, obs10, d10] = ksOfBinnedCon(chl_Out,pBin10,10);               % 10 dbar / Chlorophyll a (Regular Method)
[ksHp5, obsHp5, dHp5] = ksOfBinnedCon(chlOut_hplc,pb5_hplc,5);       % 5 dbar / HPLC Chlorophyll a
[ksHp10, obsHp10, dHp10] = ksOfBinnedCon(chlOut_hplc,pb10_hplc,10);  % 10 dbar / HPLC Chlorophyll a

[ksAtp5, obsAtp5, dAtp5] = ksOfBinnedCon(atpOut,pb5_atp,5);          % 5 dbar / ATP
[ksAtp10, obsAtp10, dAtp10] = ksOfBinnedCon(atpOut,pb10_atp,10);     % 10 dbar / ATP

[ksCmo5, obsCmo5, dCmo5] = ksOfBinnedCon(cmoOut,pb5_cmo,5);          % 5 dbar / HPLC Monovinyl Chlorophyll a
[ksCmo10, obsCmo10, dCmo10] = ksOfBinnedCon(cmoOut,pb10_cmo,10);     % 10 dbar / HPLC Monovinyl Chlorophyll a

[ksCdi5, obsCdi5, dCdi5] = ksOfBinnedCon(cdiOut,pb5_cdi,5);          % 5 dbar / HPLC Divinyl Chlorophyll a
[ksCdi10, obsCdi10, dCdi10] = ksOfBinnedCon(cdiOut,pb10_cdi,10);     % 10 dbar / HPLC Divinyl Chlorophyll a

[ksParc5, obsParc5, dParc5] = ksOfBinnedCon(parcOut,pb5_parc,5);     % 5 dbar / Particulate Carbon
[ksParc10, obsParc10, dParc10] = ksOfBinnedCon(parcOut,pb10_parc,10);% 10 dbar / Particulate Carbon

[ksNit5, obsNit5, dNit5] = ksOfBinnedCon(nitOut,pb5_nit,5);          % 5 dbar / Particulate Nitrogen
[ksNit10, obsNit10, dNit10] = ksOfBinnedCon(nitOut,pb10_nit,10);     % 10 dbar / Particulate Nitrogen

[ksPho5, obsPho5, dPho5] = ksOfBinnedCon(phoOut,pb5_pho,5);          % 5 dbar / Particulate Phosphorus
[ksPho10, obsPho10, dPho10] = ksOfBinnedCon(phoOut,pb10_pho,10);     % 10 dbar / Particulate Phosphorus

[ksPho5_11, obsPho5_11, dPho5_11] = ksOfBinnedCon(phoOut11,pb5_pho11,5);          % 5 dbar / Particulate Phosphorus 11-21
[ksPho10_11, obsPho10_11, dPho10_11] = ksOfBinnedCon(phoOut11,pb10_pho11,10);     % 10 dbar / Particulate Phosphorus 11-21

[ksHet5, obsHet5, dHet5] = ksOfBinnedCon(hetOut,pb5_het,5);          % 5 dbar / Heterotrophic Bacteria
[ksHet10, obsHet10, dHet10] = ksOfBinnedCon(hetOut,pb10_het,10);     % 10 dbar / Heterotrophic Bacteria

[ksPro5, obsPro5, dPro5] = ksOfBinnedCon(proOut,pb5_pro,5);          % 5 dbar / Prochlorococcus
[ksPro10, obsPro10, dPro10] = ksOfBinnedCon(proOut,pb10_pro,10);     % 10 dbar / Prochlorococcus

[ksSyn5, obsSyn5, dSyn5] = ksOfBinnedCon(synOut,pb5_syn,5);          % 5 dbar / Synechococcus
[ksSyn10, obsSyn10, dSyn10] = ksOfBinnedCon(synOut,pb10_syn,10);     % 10 dbar / Synechococcus

[ksPic5, obsPic5, dPic5] = ksOfBinnedCon(picOut,pb5_pic,5);          % 5 dbar / Pico-Eukaryotes
[ksPic10, obsPic10, dPic10] = ksOfBinnedCon(picOut,pb10_pic,10);     % 10 dbar / Pico-Eukaryotes

%% Chlorophyll a (regular method): Eulerian 5 dbar

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

%% HPLC Chlorophyll a: Eulerian 5 dbar

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

%% Chlorophyll a: Eulerian 10 dbar

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

%% HPLC Chlorophyll a: Eulerian 10 dbar

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

%% ATP: Eulerian 5 dbar 

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

%% ATP: Eulerian 10 dbar

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

%% HPLC Monovinyl Chlorophyll a: Eulerian 5 dbar 

ax7 = figure;
subplot(1,2,1)
barh(obsCmo5);
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
plot(ksCmo5(1,:),dCmo5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksCmo5(2,:),dCmo5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksCmo5(3,:),dCmo5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksCmo5(4,:),dCmo5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('HPLC Monovinyl Chlorophyll a: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax7,'figures/ks_CmoEulerian5db.png');
clear ax7;

%% HPLC Monovinyl Chlorophyll a: Eulerian 10 dbar

ax8 = figure;
subplot(1,2,1)
barh(obsCmo10);
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
plot(ksCmo10(1,:),dCmo10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksCmo10(2,:),dCmo10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksCmo10(3,:),dCmo10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksCmo10(4,:),dCmo10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('HPLC Monovinyl Chlorophyll a: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax8,'figures/ks_CmoEulerian10db.png');
clear ax8;

%% HPLC Divinyl Chlorophyll a: Eulerian 5 dbar 

ax9 = figure;
subplot(1,2,1)
barh(obsCdi5);
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
plot(ksCdi5(1,:),dCdi5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksCdi5(2,:),dCdi5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksCdi5(3,:),dCdi5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksCdi5(4,:),dCdi5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('HPLC Divinyl Chlorophyll a: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax9,'figures/ks_CdiEulerian5db.png');
clear ax9;

%% HPLC Divinyl Chlorophyll a: Eulerian 10 dbar

ax10 = figure;
subplot(1,2,1)
barh(obsCdi10);
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
plot(ksCdi10(1,:),dCdi10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksCdi10(2,:),dCdi10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksCdi10(3,:),dCdi10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksCdi10(4,:),dCdi10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('HPLC Divinyl Chlorophyll a: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax10,'figures/ks_CdiEulerian10db.png');
clear ax10;

%% Particulate Carbon: Eulerian 5 dbar 

ax11 = figure;
subplot(1,2,1)
barh(obsParc5);
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
plot(ksParc5(1,:),dParc5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksParc5(2,:),dParc5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksParc5(3,:),dParc5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksParc5(4,:),dParc5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Particulate Carbon: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax11,'figures/ks_ParcEulerian5db.png');
clear ax11;

%% Particulate Carbon: Eulerian 10 dbar

ax12 = figure;
subplot(1,2,1)
barh(obsParc10);
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
plot(ksParc10(1,:),dParc10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksParc10(2,:),dParc10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksParc10(3,:),dParc10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksParc10(4,:),dParc10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Particulate Carbon: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax12,'figures/ks_ParcEulerian10db.png');
clear ax12;

%% Particulate Nitrogen: Eulerian 5 dbar 

ax13 = figure;
subplot(1,2,1)
barh(obsNit5);
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
plot(ksNit5(1,:),dNit5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksNit5(2,:),dNit5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksNit5(3,:),dNit5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksNit5(4,:),dNit5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Particulate Nitrogen: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax13,'figures/ks_NitEulerian5db.png');
clear ax13;

%% Particulate Nitrogen: Eulerian 10 dbar

ax14 = figure;
subplot(1,2,1)
barh(obsNit10);
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
plot(ksNit10(1,:),dNit10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksNit10(2,:),dNit10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksNit10(3,:),dNit10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksNit10(4,:),dNit10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Particulate Nitrogen: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax14,'figures/ks_NitEulerian10db.png');
clear ax14;

%% Particulate Phosphorus: Eulerian 5 dbar 

ax15 = figure;
subplot(1,2,1)
barh(obsPho5);
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
plot(ksPho5(1,:),dPho5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPho5(2,:),dPho5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPho5(3,:),dPho5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPho5(4,:),dPho5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Particulate Phosphorus: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax15,'figures/ks_PhoEulerian5db.png');
clear ax15;

%% Particulate Phosphorus: Eulerian 10 dbar

ax16 = figure;
subplot(1,2,1)
barh(obsPho10);
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
plot(ksPho10(1,:),dPho10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPho10(2,:),dPho10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPho10(3,:),dPho10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPho10(4,:),dPho10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Particulate Phosphorus: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax16,'figures/ks_PhoEulerian10db.png');
clear ax16;

%% Heterotrophic Bacteria: Eulerian 5 dbar 

ax17 = figure;
subplot(1,2,1)
barh(obsHet5);
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
plot(ksHet5(1,:),dHet5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksHet5(2,:),dHet5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksHet5(3,:),dHet5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksHet5(4,:),dHet5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Heterotrophic Bacteria: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax17,'figures/ks_HetEulerian5db.png');
clear ax17;

%% Particulate Phosphorus: Eulerian 10 dbar

ax18 = figure;
subplot(1,2,1)
barh(obsHet10);
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
plot(ksHet10(1,:),dHet10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksHet10(2,:),dHet10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksHet10(3,:),dHet10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksHet10(4,:),dHet10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Heterotrophic Bacteria: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax18,'figures/ks_HetEulerian10db.png');
clear ax18;

%% Prochlorococcus: Eulerian 5 dbar 

ax19 = figure;
subplot(1,2,1)
barh(obsPro5);
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
plot(ksPro5(1,:),dPro5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPro5(2,:),dPro5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPro5(3,:),dPro5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPro5(4,:),dPro5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Prochlorococcus: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax19,'figures/ks_ProEulerian5db.png');
clear ax19;

%% Prochlorococcus: Eulerian 10 dbar

ax20 = figure;
subplot(1,2,1)
barh(obsPro10);
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
plot(ksPro10(1,:),dPro10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPro10(2,:),dPro10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPro10(3,:),dPro10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPro10(4,:),dPro10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Prochlorococcus: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax20,'figures/ks_ProEulerian10db.png');
clear ax20;

%% Synechococcus: Eulerian 5 dbar 

ax21 = figure;
subplot(1,2,1)
barh(obsSyn5);
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
plot(ksSyn5(1,:),dSyn5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksSyn5(2,:),dSyn5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksSyn5(3,:),dSyn5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksSyn5(4,:),dSyn5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Synechococcus: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax21,'figures/ks_SynEulerian5db.png');
clear ax21;

%% Synechococcus: Eulerian 10 dbar

ax22 = figure;
subplot(1,2,1)
barh(obsSyn10);
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
plot(ksSyn10(1,:),dSyn10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksSyn10(2,:),dSyn10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksSyn10(3,:),dSyn10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksSyn10(4,:),dSyn10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Synechococcus: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax22,'figures/ks_SynEulerian10db.png');
clear ax22;

%% Pico-eukaryotes: Eulerian 5 dbar 

ax23 = figure;
subplot(1,2,1)
barh(obsPic5);
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
plot(ksPic5(1,:),dPic5,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPic5(2,:),dPic5,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPic5(3,:),dPic5,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPic5(4,:),dPic5,'r.--','DisplayName','Gamma','MarkerSize',6);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Pico-Eukaryotes: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax23,'figures/ks_PicEulerian5db.png');
clear ax23;

%% Pico Eukaryotes: Eulerian 10 dbar

ax24 = figure;
subplot(1,2,1)
barh(obsPic10);
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
plot(ksPic10(1,:),dPic10,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPic10(2,:),dPic10,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPic10(3,:),dPic10,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPic10(4,:),dPic10,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('Kolmogorov Smirnov Test');

sgtitle('Pico-Eukaryotes: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax24,'figures/ks_PicEulerian10db.png');
clear ax24;

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