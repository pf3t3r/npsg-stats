clear; clc; close all; addpath("baroneRoutines\");

%% Load MLD and DCM
pMaxMld = load('mldVals.mat').maxMld;
dcm = load("dcm.mat").dcm;  % pDcm and sigmaDcm (all casts, crn 1 - 329)

%% Load parameters to test
% Bottle ID, pressure, and concentration

% Chlorophyll a
chlId = num2str(importdata('data/L2/hplcChla_88-21_200.txt').data(:,1));
chlP = importdata('data/L2/hplcChla_88-21_200.txt').data(:,4);
chl = importdata('data/L2/hplcChla_88-21_200.txt').data(:,5);

% Divinyl Chlorophyll a
divId = num2str(importdata('data/L2/diviChla_88-21_200.txt').data(:,1));
divP = importdata('data/L2/diviChla_88-21_200.txt').data(:,4);
div = importdata('data/L2/diviChla_88-21_200.txt').data(:,5);

% Prochlorococcus: 05-21
proId = num2str(importdata('data/L2/pro_05-21_200.txt').data(:,1));
proP = importdata('data/L2/pro_05-21_200.txt').data(:,4);
pro = importdata('data/L2/pro_05-21_200.txt').data(:,5);

% Prochlorococcus: 90-05
proId9 = num2str(importdata('data/L2/pro_90-05_200.txt').data(:,1));
proP9 = importdata('data/L2/pro_90-05_200.txt').data(:,4);
pro9 = importdata('data/L2/pro_90-05_200.txt').data(:,5);

% Oxygen: 88-21
oxyId = num2str(importdata('data/L2/oxy_88-21_200.txt').data(:,1));
oxyP = importdata('data/L2/oxy_88-21_200.txt').data(:,4);
oxy = importdata('data/L2/oxy_88-21_200.txt').data(:,5);

% Dissolved Inorganic Carbon: 88-21
dicId = num2str(importdata('data/L2/dic_88-21_200.txt').data(:,1));
dicP = importdata('data/L2/dic_88-21_200.txt').data(:,4);
dic = importdata('data/L2/dic_88-21_200.txt').data(:,5);

% pH: 92-21
phId = num2str(importdata('data/L2/pH_92-21_200.txt').data(:,1));
phP = importdata('data/L2/pH_92-21_200.txt').data(:,4);
pH = importdata('data/L2/pH_92-21_200.txt').data(:,5);

% Alkalinity: 89-21
alkId = num2str(importdata('data/L2/alk_89-21_200.txt').data(:,1));
alkP = importdata('data/L2/alk_89-21_200.txt').data(:,4);
alk = importdata('data/L2/alk_89-21_200.txt').data(:,5);

% Phosphate: 88-21
phoId = num2str(importdata('data/L2/pho_88-21_200.txt').data(:,1));
phoP = importdata('data/L2/pho_88-21_200.txt').data(:,4);
pho = importdata('data/L2/pho_88-21_200.txt').data(:,5);

% Nitrate + Nitrite: 88-21
nit2Id = num2str(importdata('data/L2/nit2_88-21_200.txt').data(:,1));
nit2P = importdata('data/L2/nit2_88-21_200.txt').data(:,4);
nit2 = importdata('data/L2/nit2_88-21_200.txt').data(:,5);

% Nitrite: 89-94
nitId = num2str(importdata('data/L2/nit_89-94_200.txt').data(:,1));
nitP = importdata('data/L2/nit_89-94_200.txt').data(:,4);
nit = importdata('data/L2/nit_89-94_200.txt').data(:,5);

%% Extract sub-ML bottle ID, pressure, bottle concentrations

% extractSMLC = 'extract sub-ML concentration'
[idSubmlChl,pSubmlChl,submlChl] = extractSMLC(chlId,chlP,chl,pMaxMld);
[idSubmlDiv,pSubmlDiv,submlDiv] = extractSMLC(divId,divP,div,pMaxMld);
[idSubmlPro,pSubmlPro,submlPro] = extractSMLC(proId,proP,pro,pMaxMld);
[idSubmlPro9,pSubmlPro9,submlPro9] = extractSMLC(proId9,proP9,pro9,pMaxMld);
[idSubmlOxy,pSubmlOxy,submlOxy] = extractSMLC(oxyId,oxyP,oxy,pMaxMld);
[idSubmlDic,pSubmlDic,submlDic] = extractSMLC(dicId,dicP,dic,pMaxMld);
[idSubmlPh,pSubmlPh,submlPh] = extractSMLC(phId,phP,pH,pMaxMld);
[idSubmlAlk,pSubmlAlk,submlAlk] = extractSMLC(alkId,alkP,alk,pMaxMld);
[idSubmlPho,pSubmlPho,submlPho] = extractSMLC(phoId,phoP,pho,pMaxMld);
[idSubmlNit2,pSubmlNit2,submlNit2] = extractSMLC(nit2Id,nit2P,nit2,pMaxMld);
% [idSubmlNit,pSubmlNit,submlNit] = extractSMLC(nitId,nitP,nit,pMaxMld);

%% KS Test p-values, skewness, kurtosis: calculate

[prChl,ksChl,obsChl,skChl,kuChl,~,~,~] = ksOfLagrangian(idSubmlChl,pSubmlChl,dcm,submlChl,81);
[prDiv,ksDiv,obsDiv,skDiv,kuDiv,~,~,~] = ksOfLagrangian(idSubmlDiv,pSubmlDiv,dcm,submlDiv,69);
[prPro,ksPro,obsPro,skPro,kuPro,~,~,~] = ksOfLagrangian(idSubmlPro,pSubmlPro,dcm,submlPro,42);
[prPro9,ksPro9,obsPro9,skPro9,kuPro9,~,~,~] = ksOfLagrangian(idSubmlPro9,pSubmlPro9,dcm,submlPro9,30);
[prOxy,ksOxy,obsOxy,skOxy,kuOxy,~,~,~] = ksOfLagrangian(idSubmlOxy,pSubmlOxy,dcm,submlOxy,122);
[prDic,ksDic,obsDic,skDic,kuDic,~,~,~] = ksOfLagrangian(idSubmlDic,pSubmlDic,dcm,submlDic,89);
[prPh,ksPh,obsPh,skPh,kuPh,~,~,~] = ksOfLagrangian(idSubmlPh,pSubmlPh,dcm,submlPh,59);
[prAlk,ksAlk,obsAlk,skAlk,kuAlk,~,~,~] = ksOfLagrangian(idSubmlAlk,pSubmlAlk,dcm,submlAlk,83);
[prPho,ksPho,obsPho,skPho,kuPho,~,~,~] = ksOfLagrangian(idSubmlPho,pSubmlPho,dcm,submlPho);
[prNit2,ksNit2,obsNit2,skNit2,kuNit2,~,~,~] = ksOfLagrangian(idSubmlNit2,pSubmlNit2,dcm,submlNit2);
% [prNit,ksNit,obsNit,skNit,kuNit,~,~,~] = ksOfLagrangian(idSubmlNit,pSubmlNit,dcm,submlNit);

%% KS Test p-values, skewness, kurtosis: visualise

ax1 = figure;
plotKs2(prChl,ksChl,obsChl,skChl,kuChl,prChl(1),prChl(end),81);
sgtitle('[chl a]');
exportgraphics(ax1,'figures/L2/ks_chla.png'); clear ax1;

ax2 = figure;
plotKs2(prDiv,ksDiv,obsDiv,skDiv,kuDiv,prDiv(1),prDiv(end),69);
sgtitle('[Divinyl Chl a]');
exportgraphics(ax2,'figures/L2/ks_div.png'); clear ax2;

ax3 = figure;
plotKs2(prPro,ksPro,obsPro,skPro,kuPro,prPro(1),prPro(end),42);
sgtitle('Prochlorococcus: 05-21');
exportgraphics(ax3,'figures/L2/ks_pro.png'); clear ax3;

ax4 = figure;
plotKs2(prPro9,ksPro9,obsPro9,skPro9,kuPro9,prPro9(1),prPro9(end),30);
sgtitle('Prochlorococcus: 90-05');
exportgraphics(ax4,'figures/L2/ks_pro9.png'); clear ax4;

ax5 = figure;
plotKs2(prOxy,ksOxy,obsOxy,skOxy,kuOxy,prOxy(1),prOxy(end),122);
sgtitle('Dissolved Oxygen: 88-21');
exportgraphics(ax5,'figures/L2/ks_oxy.png'); clear ax5;

ax6 = figure;
plotKs2(prDic,ksDic,obsDic,skDic,kuDic,prDic(1),prDic(end),89);
sgtitle('Dissolved Inorganic Carbon: 88-21');
exportgraphics(ax6,'figures/L2/ks_dic.png'); clear ax6;

ax7 = figure;
plotKs2(prPh,ksPh,obsPh,skPh,kuPh,prPh(1),prPh(end),59);
sgtitle('pH: 92-21');
exportgraphics(ax7,'figures/L2/ks_pH.png'); clear ax7;

ax8 = figure;
plotKs2(prAlk,ksAlk,obsAlk,skAlk,kuAlk,prAlk(1),prAlk(end),83);
sgtitle('Alkalinity: 89-21');
exportgraphics(ax8,'figures/L2/ks_alk.png'); clear ax8;

ax9 = figure;
plotKs2(prPho,ksPho,obsPho,skPho,kuPho,prPho(1),prPho(end));
sgtitle('Phosphate: 88-21');
exportgraphics(ax9,'figures/L2/ks_pho.png'); clear ax9;

ax10 = figure;
plotKs2(prNit2,ksNit2,obsNit2,skNit2,kuNit2,prNit2(1),prNit2(end));
sgtitle('Nitrate + Nitrite: 88-21');
exportgraphics(ax10,'figures/L2/ks_nit2.png'); clear ax10;

% ax11 = figure;
% plotKs2(prNit,ksNit,obsNit,skNit,kuNit,prNit(1),prNit(end));
% sgtitle('Nitrite: 89-94');
% exportgraphics(ax11,'figures/L2/ks_nit.png'); clear ax11;