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
[pb5_pho11,pb10_pho11,phoOut11,n5_pho11,n10_pho11] = cleanAndBin(p_pho11,pho11,id_pho11);       % Particulate Phosphorus (11-21)
[pb5_het,pb10_het,hetOut,n5_het,n10_het] = cleanAndBin(p_het,het,id_het);                       % Heterotrophic Bacteria
[pb5_pro,pb10_pro,proOut,n5_pro,n10_pro] = cleanAndBin(p_pro,pro,id_pro);                       % Prochlorococcus
[pb5_syn,pb10_syn,synOut,n5_syn,n10_syn] = cleanAndBin(p_syn,syn,id_syn);                       % Synechococcus
[pb5_pic,pb10_pic,picOut,n5_pic,n10_pic] = cleanAndBin(p_pic,pic,id_pic);                       % Pico-Eukaryotes

%% Apply KS Test to chl-a across all pressures

[ks5, obs5, d5, sk5, ku5, sd5, c95_5, mu5] = ksOfBinnedCon(chl_Out,pBin5,5);                    % 5 dbar / Chlorophyll a (Regular Method)
[ks10, obs10, d10, sk10, ku10, sd10, c95_10, mu10] = ksOfBinnedCon(chl_Out,pBin10,10);               % 10 dbar / Chlorophyll a (Regular Method)
[ksHp5, obsHp5, dHp5, skHp5, kuHp5, sdHp5, c95hp5, muHp5] = ksOfBinnedCon(chlOut_hplc,pb5_hplc,5);       % 5 dbar / HPLC Chlorophyll a
[ksHp10, obsHp10, dHp10, skHp10, kuHp10, sdHp10, c95hp10, muHp10] = ksOfBinnedCon(chlOut_hplc,pb10_hplc,10);  % 10 dbar / HPLC Chlorophyll a

[ksAtp5, obsAtp5, dAtp5, skAtp5, kuAtp5, sdAtp5, c95atp5, muAtp5] = ksOfBinnedCon(atpOut,pb5_atp,5);          % 5 dbar / ATP
[ksAtp10, obsAtp10, dAtp10, skAtp10, kuAtp10, sdAtp10, c95atp10, muAtp10] = ksOfBinnedCon(atpOut,pb10_atp,10);     % 10 dbar / ATP

[ksCmo5, obsCmo5, dCmo5, skCmo5, kuCmo5, sdCmo5, c95cmo5, muCmo5] = ksOfBinnedCon(cmoOut,pb5_cmo,5);          % 5 dbar / HPLC Monovinyl Chlorophyll a
[ksCmo10, obsCmo10, dCmo10, skCmo10, kuCmo10, sdCmo10, c95cmo10, muCmo10] = ksOfBinnedCon(cmoOut,pb10_cmo,10);     % 10 dbar / HPLC Monovinyl Chlorophyll a

[ksCdi5, obsCdi5, dCdi5, skCdi5, kuCdi5, sdCdi5, c95cdi5, muCdi5] = ksOfBinnedCon(cdiOut,pb5_cdi,5);          % 5 dbar / HPLC Divinyl Chlorophyll a
[ksCdi10, obsCdi10, dCdi10, skCdi10, kuCdi10, sdCdi10, c95cdi10, muCdi10] = ksOfBinnedCon(cdiOut,pb10_cdi,10);     % 10 dbar / HPLC Divinyl Chlorophyll a

[ksParc5, obsParc5, dParc5, skParc5, kuParc5, sdParc5, c95parc5, muParc5] = ksOfBinnedCon(parcOut,pb5_parc,5);     % 5 dbar / Particulate Carbon
[ksParc10, obsParc10, dParc10, skParc10, kuParc10, sdParc10, c95parc10, muParc10] = ksOfBinnedCon(parcOut,pb10_parc,10);% 10 dbar / Particulate Carbon

[ksNit5, obsNit5, dNit5, skNit5, kuNit5, sdNit5, c95nit5, muNit5] = ksOfBinnedCon(nitOut,pb5_nit,5);          % 5 dbar / Particulate Nitrogen
[ksNit10, obsNit10, dNit10, skNit10, kuNit10, sdNit10, c95nit10, muNit10] = ksOfBinnedCon(nitOut,pb10_nit,10);     % 10 dbar / Particulate Nitrogen

[ksPho5, obsPho5, dPho5, skPho5, kuPho5, sdPho5, c95pho5, muPho5] = ksOfBinnedCon(phoOut,pb5_pho,5);          % 5 dbar / Particulate Phosphorus
[ksPho10, obsPho10, dPho10, skPho10, kuPho10, sdPho10, c95pho10, muPho10] = ksOfBinnedCon(phoOut,pb10_pho,10);     % 10 dbar / Particulate Phosphorus

[ksPho5_11, obsPho5_11, dPho5_11, skPho5_11, kuPho5_11, sdPho5_11, c95pho5_11, muPho5_11] = ksOfBinnedCon(phoOut11,pb5_pho11,5,75);          % 5 dbar / Particulate Phosphorus 11-21
[ksPho10_11, obsPho10_11, dPho10_11, skPho10_11, kuPho10_11, sdPho10_11, c95pho10_11, muPho10_11] = ksOfBinnedCon(phoOut11,pb10_pho11,10,75);     % 10 dbar / Particulate Phosphorus 11-21

[ksHet5, obsHet5, dHet5, skHet5, kuHet5, sdHet5, c95het5, muHet5] = ksOfBinnedCon(hetOut,pb5_het,5);          % 5 dbar / Heterotrophic Bacteria
[ksHet10, obsHet10, dHet10, skHet10, kuHet10, sdHet10, c95het10, muHet10] = ksOfBinnedCon(hetOut,pb10_het,10);     % 10 dbar / Heterotrophic Bacteria

[ksPro5, obsPro5, dPro5, skPro5, kuPro5, sdPro5, c95Pro5, muPro5] = ksOfBinnedCon(proOut,pb5_pro,5);          % 5 dbar / Prochlorococcus
[ksPro10, obsPro10, dPro10, skPro10, kuPro10, sdPro10, c95Pro10, muPro10] = ksOfBinnedCon(proOut,pb10_pro,10);     % 10 dbar / Prochlorococcus

[ksSyn5, obsSyn5, dSyn5, skSyn5, kuSyn5, sdSyn5, c95syn5, muSyn5] = ksOfBinnedCon(synOut,pb5_syn,5);          % 5 dbar / Synechococcus
[ksSyn10, obsSyn10, dSyn10, skSyn10, kuSyn10, sdSyn10, c95syn10, muSyn10] = ksOfBinnedCon(synOut,pb10_syn,10);     % 10 dbar / Synechococcus

[ksPic5, obsPic5, dPic5, skPic5, kuPic5, sdPic5, c95pic5, muPic5] = ksOfBinnedCon(picOut,pb5_pic,5);          % 5 dbar / Pico-Eukaryotes
[ksPic10, obsPic10, dPic10, skPic10, kuPic10, sdPic10, c95pic10, muPic10] = ksOfBinnedCon(picOut,pb10_pic,10);     % 10 dbar / Pico-Eukaryotes

%% Chlorophyll a (regular method): Eulerian 5 dbar

ax1 = figure;
plotKs(d5,ks5,obs5,sk5,ku5,0,40,true);
sgtitle('Regular Method: [chl a] (Eulerian Bottle, 5 dbar bin)');
exportgraphics(ax1,'figures/ks_bottleEulerian5db.png'); clear ax1;

%% HPLC Chlorophyll a: Eulerian 5 dbar

ax2 = figure;
plotKs(dHp5,ksHp5,obsHp5,skHp5,kuHp5,0,40,true);
sgtitle('HPLC Method: [chl a] (Eulerian Bottle, 5 dbar bin)');
exportgraphics(ax2,'figures/ks_HplcBottleEulerian5db.png'); clear ax2;

%% Chlorophyll a: Eulerian 10 dbar

ax3 = figure;
plotKs(d10,ks10,obs10,sk10,ku10,0.5,20.5,true);
sgtitle('Regular Method: [chl a] (Eulerian Bottle, 10 dbar bin)');
exportgraphics(ax3,'figures/ks_bottleEulerian10db.png'); clear ax3;

%% HPLC Chlorophyll a: Eulerian 10 dbar

ax4 = figure;
plotKs(dHp10,ksHp10,obsHp10,skHp10,kuHp10,0.5,20.5,true);
sgtitle('HPLC Method: [chl a] (Eulerian Bottle, 10 dbar bin)');
exportgraphics(ax4,'figures/ks_HplcBottleEulerian10db.png'); clear ax4;

%% ATP: Eulerian 5 dbar 

ax5 = figure;
plotKs(dAtp5,ksAtp5,obsAtp5,skAtp5,kuAtp5,0,40,true);
sgtitle('ATP: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax5,'figures/ks_AtpEulerian5db.png'); clear ax5;

%% ATP: Eulerian 10 dbar

ax6 = figure;
plotKs(dAtp10,ksAtp10,obsAtp10,skAtp10,kuAtp10,0.5,20.5,true);
sgtitle('ATP: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax6,'figures/ks_AtpEulerian10db.png'); clear ax6;

%% HPLC Monovinyl Chlorophyll a: Eulerian 5 dbar 

ax7 = figure;
plotKs(dCmo5,ksCmo5,obsCmo5,skCmo5,kuCmo5,0,40,true);
sgtitle('HPLC Monovinyl Chlorophyll a: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax7,'figures/ks_CmoEulerian5db.png'); clear ax7;

%% HPLC Monovinyl Chlorophyll a: Eulerian 10 dbar

ax8 = figure;
plotKs(dCmo10,ksCmo10,obsCmo10,skCmo10,kuCmo10,0.5,20.5,true);
sgtitle('HPLC Monovinyl Chlorophyll a: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax8,'figures/ks_CmoEulerian10db.png'); clear ax8;

%% HPLC Divinyl Chlorophyll a: Eulerian 5 dbar 

ax9 = figure;
plotKs(dCdi5,ksCdi5,obsCdi5,skCdi5,kuCdi5,0,40,true);
sgtitle('HPLC Divinyl Chlorophyll a: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax9,'figures/ks_CdiEulerian5db.png'); clear ax9;

%% HPLC Divinyl Chlorophyll a: Eulerian 10 dbar

ax10 = figure;
plotKs(dCdi10,ksCdi10,obsCdi10,skCdi10,kuCdi10,0.5,20.5,true);
sgtitle('HPLC Divinyl Chlorophyll a: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax10,'figures/ks_CdiEulerian10db.png'); clear ax10;

%% Particulate Carbon: Eulerian 5 dbar 

ax11 = figure;
plotKs(dParc5,ksParc5,obsParc5,skParc5,kuParc5,0,40,true);
sgtitle('Particulate Carbon: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax11,'figures/ks_ParcEulerian5db.png'); clear ax11;

%% Particulate Carbon: Eulerian 10 dbar

ax12 = figure;
plotKs(dParc10,ksParc10,obsParc10,skParc10,kuParc10,0.5,20.5,true);
sgtitle('Particulate Carbon: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax12,'figures/ks_ParcEulerian10db.png'); clear ax12;

%% Particulate Nitrogen: Eulerian 5 dbar 

ax13 = figure;
plotKs(dNit5,ksNit5,obsNit5,skNit5,kuNit5,0,40,true);
sgtitle('Particulate Nitrogen: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax13,'figures/ks_NitEulerian5db.png'); clear ax13;

%% Particulate Nitrogen: Eulerian 10 dbar

ax14 = figure;
plotKs(dNit10,ksNit10,obsNit10,skNit10,kuNit10,0.5,20.5,true);
sgtitle('Particulate Nitrogen: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax14,'figures/ks_NitEulerian10db.png'); clear ax14;

%% Particulate Phosphorus: Eulerian 5 dbar 

ax15 = figure;
plotKs(dPho5,ksPho5,obsPho5,skPho5,kuPho5,0,40,true);
sgtitle('Particulate Phosphorus: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax15,'figures/ks_PhoEulerian5db.png'); clear ax15;

ax15a = figure;
plotKs(dPho5_11,ksPho5_11,obsPho5_11,skPho5_11,kuPho5_11,0,40,true,75);
sgtitle('Particulate Phosphorus: Eulerian Bottle, Threshold=75, 5 dbar bin (2011-2021)');
exportgraphics(ax15a,'figures/ks_PhoEulerian5db_11-21thr75.png'); clear ax15a;

%% Particulate Phosphorus: Eulerian 10 dbar

ax16 = figure;
plotKs(dPho10,ksPho10,obsPho10,skPho10,kuPho10,0.5,20.5,true);
sgtitle('Particulate Phosphorus: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax16,'figures/ks_PhoEulerian10db.png'); clear ax16;

ax16a = figure;
plotKs(dPho10_11,ksPho10_11,obsPho10_11,skPho10_11,kuPho10_11,0.5,20.5,true,75);
sgtitle('Particulate Phosphorus: Eulerian Bottle, Threshold=75, 10 dbar bin (2011-2021)');
exportgraphics(ax16a,'figures/ks_PhoEulerian10db_11-21thr75.png'); clear ax16a;

%% plotting xx example

axTest = figure;
subplot(1,4,3)
plot(skPho10,dPho10,'Color','#1f78b4');
tAx = gca;
set(tAx,'YDir','reverse');
tAx2 = axes('position', get(tAx, 'position')); % Create a new axes in the same position as the first one, overlaid on top
plot(kuPho10,dPho10,'Color','#33a02c');
set(tAx2, 'ylim', get(tAx, 'ylim'), 'color', 'none'); % Set y limits same as original axes, and make background transparent
set(tAx2,'XAxisLocation','top');
ylabel(tAx, 'Contour y-values');
xlabel(tAx,'skew','Color','#1f78b4');
xlabel(tAx2,'kurt','Color','#33a02c');
set(tAx2,'YTickLabel',[]);
set(tAx,'XColor','#1f78b4');
set(tAx2,'XColor','#33a02c');
set(tAx2,'YDir','reverse');

% clear axTest tAx tAx2;

%% Heterotrophic Bacteria: Eulerian 5 dbar 

ax17 = figure;
plotKs(dHet5,ksHet5,obsHet5,skHet5,kuHet5,0,40,true);
sgtitle('Heterotrophic Bacteria: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax17,'figures/ks_HetEulerian5db.png'); clear ax17;

%% Het Bac: Eulerian 10 dbar

ax18 = figure;
plotKs(dHet10,ksHet10,obsHet10,skHet10,kuHet10,0.5,20.5,true);
sgtitle('Heterotrophic Bacteria: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax18,'figures/ks_HetEulerian10db.png'); clear ax18;

%% Prochlorococcus: Eulerian 5 dbar 

ax19 = figure;
plotKs(dPro5,ksPro5,obsPro5,skPro5,kuPro5,0,40,true);
sgtitle('Prochlorococcus: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax19,'figures/ks_ProEulerian5db.png'); clear ax19;

%% Prochlorococcus: Eulerian 10 dbar

ax20 = figure;
plotKs(dPro10,ksPro10,obsPro10,skPro10,kuPro10,0.5,20.5,true);
sgtitle('Prochlorococcus: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax20,'figures/ks_ProEulerian10db.png'); clear ax20;

%% Synechococcus: Eulerian 5 dbar 

ax21 = figure;
plotKs(dSyn5,ksSyn5,obsSyn5,skSyn5,kuSyn5,0,40,true);
sgtitle('Synechococcus: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax21,'figures/ks_SynEulerian5db.png'); clear ax21;

%% Synechococcus: Eulerian 10 dbar

ax22 = figure;
plotKs(dSyn10,ksSyn10,obsSyn10,skSyn10,kuSyn10,0.5,20.5,true);
sgtitle('Synechococcus: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax22,'figures/ks_SynEulerian10db.png'); clear ax22;

%% Pico-eukaryotes: Eulerian 5 dbar 

ax23 = figure;
plotKs(dPic5,ksPic5,obsPic5,skPic5,kuPic5,0,40,true);
sgtitle('Pico-Eukaryotes: Eulerian Bottle, 5 dbar bin');
exportgraphics(ax23,'figures/ks_PicEulerian5db.png'); clear ax23;

%% Pico Eukaryotes: Eulerian 10 dbar

ax24 = figure;
plotKs(dPic10,ksPic10,obsPic10,skPic10,kuPic10,0.5,20.5,true);
sgtitle('Pico-Eukaryotes: Eulerian Bottle, 10 dbar bin');
exportgraphics(ax24,'figures/ks_PicEulerian10db.png'); clear ax24;

%% Visualise STD Fraction

norm = 5; logn = 6;

xax = 1:1:20;
ax25 = figure;
subplot(2,1,1)
scatter(sd10(:,norm),d10,'DisplayName','[chl a] (fluo.)','Marker','+');
hold on
scatter(sdHp10(:,norm),dHp10,'DisplayName','[chl a] (HPLC)','Marker','+');
scatter(sdAtp10(:,norm),dAtp10,'DisplayName','ATP','Marker','+');
scatter(sdCmo10(:,norm),dCmo10,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
scatter(sdCdi10(:,norm),dCdi10,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
scatter(sdParc10(:,norm),dParc10,'DisplayName','Particulate Carbon','Marker','+');
scatter(sdNit10(:,norm),dNit10,'DisplayName','Particulate Nitrogen','Marker','+');
scatter(sdPho10(:,norm),dPho10,'DisplayName','Particulate Phosphorus');
scatter(sdHet10(:,norm),dHet10,'DisplayName','Heterotrophic Bacteria ');
scatter(sdPro10(:,norm),dPro10,'DisplayName','Prochlorococcus');
scatter(sdSyn10(:,norm),dSyn10,'DisplayName','Synechococcus');
hold off
xlim([0.995 1.00]);
legend('Location','bestoutside'); set(gca,'YDir','reverse');
h = title('$\frac{\sigma_{mle}}{\sigma_{data}}$ normal','Interpreter','latex');
h.FontSize=30;

% tmpY = 5:10:200;
subplot(2,1,2)
scatter(sd10(:,logn),d10,'DisplayName','[chl a] (fluo.)','Marker','+');
hold on
scatter(sdHp10(:,logn),dHp10,'DisplayName','[chl a] (HPLC)','Marker','+');
scatter(sdAtp10(:,logn),dAtp10,'DisplayName','ATP','Marker','+');
scatter(sdCmo10(:,logn),dCmo10,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
scatter(sdCdi10(:,logn),dCdi10,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
scatter(sdParc10(:,logn),dParc10,'DisplayName','Particulate Carbon','Marker','+');
scatter(sdNit10(:,logn),dNit10,'DisplayName','Particulate Nitrogen','Marker','+');
scatter(sdPho10(:,logn),dPho10,'DisplayName','Particulate Phosphorus');
scatter(sdHet10(:,logn),dHet10,'DisplayName','Heterotrophic Bacteria ');
scatter(sdPro10(:,logn),dPro10,'DisplayName','Prochlorococcus');
scatter(sdSyn10(:,logn),dSyn10,'DisplayName','Synechococcus');
hold off
legend('Location','bestoutside'); set(gca,'YDir','reverse');
h = title('$\frac{\sigma_{mle}}{\sigma_{data}}$ lognormal','Interpreter','latex');
h.FontSize=30;

exportgraphics(ax25,'figures/stdComp.png'); clear ax25;

%% Visualise MU Fraction

norm = 5; logn = 6;

xax = 1:1:20;
ax26 = figure;
subplot(2,1,1)
scatter(mu10(:,norm),d10,'DisplayName','[chl a] (fluo.)','Marker','+');
hold on
scatter(muHp10(:,norm),dHp10,'DisplayName','[chl a] (HPLC)','Marker','+');
scatter(muAtp10(:,norm),dAtp10,'DisplayName','ATP','Marker','+');
scatter(muCmo10(:,norm),dCmo10,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
scatter(muCdi10(:,norm),dCdi10,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
scatter(muParc10(:,norm),dParc10,'DisplayName','Particulate Carbon','Marker','+');
scatter(muNit10(:,norm),dNit10,'DisplayName','Particulate Nitrogen','Marker','+');
scatter(muPho10(:,norm),dPho10,'DisplayName','Particulate Phosphorus');
scatter(muHet10(:,norm),dHet10,'DisplayName','Heterotrophic Bacteria ');
scatter(muPro10(:,norm),dPro10,'DisplayName','Prochlorococcus');
scatter(muSyn10(:,norm),dSyn10,'DisplayName','Synechococcus');
hold off
xlim([0.995 1.00]);
legend('Location','bestoutside'); set(gca,'YDir','reverse');
h = title('$\frac{\mu_{mle}}{\mu_{data}}$ normal','Interpreter','latex');
h.FontSize=30;

% tmpY = 5:10:200;
subplot(2,1,2)
scatter(mu10(:,logn),d10,'DisplayName','[chl a] (fluo.)','Marker','+');
hold on
scatter(muHp10(:,logn),dHp10,'DisplayName','[chl a] (HPLC)','Marker','+');
scatter(muAtp10(:,logn),dAtp10,'DisplayName','ATP','Marker','+');
scatter(muCmo10(:,logn),dCmo10,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
scatter(muCdi10(:,logn),dCdi10,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
scatter(muParc10(:,logn),dParc10,'DisplayName','Particulate Carbon','Marker','+');
scatter(muNit10(:,logn),dNit10,'DisplayName','Particulate Nitrogen','Marker','+');
scatter(muPho10(:,logn),dPho10,'DisplayName','Particulate Phosphorus');
scatter(muHet10(:,logn),dHet10,'DisplayName','Heterotrophic Bacteria ');
scatter(muPro10(:,logn),dPro10,'DisplayName','Prochlorococcus');
scatter(muSyn10(:,logn),dSyn10,'DisplayName','Synechococcus');
hold off
legend('Location','bestoutside'); set(gca,'YDir','reverse');
h = title('$\frac{\mu_{mle}}{\mu_{data}}$ lognormal','Interpreter','latex');
h.FontSize=30;

exportgraphics(ax26,'figures/meanComp.png'); clear ax26;


%% Show STD versus each other
% 
% figure
% scatter(sd10(:,1),sd10(:,3),'+','DisplayName','[chl a] (fluo.)');
% hold on
% scatter(sdHp10(:,1),sdHp10(:,3),'+','DisplayName','[chl a] (HPLC)');
% scatter(sdAtp10(:,1),sdAtp10(:,3),'+','DisplayName','ATP');
% scatter(sdCmo10(:,1),sdCmo10(:,3),'+','DisplayName','HPLC Monovinyl Chl');
% scatter(sdCdi10(:,1),sdCdi10(:,3),'+','DisplayName','HPLC Divinyl Chl');
% scatter(sdParc10(:,1),sdParc10(:,3),'+','DisplayName','Particulate Carbon');
% scatter(sdNit10(:,1),sdNit10(:,3),'+','DisplayName','Particulate Nitrogen');
% scatter(sdPho10(:,1),sdPho10(:,3),'o','DisplayName','Particulate Phosphorus');
% scatter(sdHet10(:,1),sdHet10(:,3),'o','DisplayName','Heterotrophic Bacteria');
% scatter(sdPro10(:,1),sdPro10(:,3),'o','DisplayName','Prochlorococcus');
% scatter(sdSyn10(:,1),sdSyn10(:,3),'o','DisplayName','Synechococcus');
% plot(1:1:100,1:1:100,':','DisplayName','slope = 1','Color',[0.6 0.6 0.6]);
% hold off
% legend('Location','bestoutside');
% xlabel('STD (MLE)'); ylabel('STD (data)');


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