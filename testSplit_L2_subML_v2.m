clear; clc; close all; addpath("baroneRoutines\");

%% Load MLD and DCM
pMaxMld = load('mldVals.mat').maxMld;
dcm = load("dcm.mat").dcm;  % pDcm and sigmaDcm (all casts, crn 1 - 329)

%% Load parameters to test
% Bottle ID, pressure, and concentration

% Chlorophyll a
chlId = num2str(importdata('data/hplcChla_88-21_200.txt').data(:,1));
chlP = importdata('data/hplcChla_88-21_200.txt').data(:,4);
chl = importdata('data/hplcChla_88-21_200.txt').data(:,5);

% Divinyl Chlorophyll a
divId = num2str(importdata('data/diviChla_88-21_200.txt').data(:,1));
divP = importdata('data/diviChla_88-21_200.txt').data(:,4);
div = importdata('data/diviChla_88-21_200.txt').data(:,5);

% Prochlorococcus: 05-21
proId = num2str(importdata('data/pro_05-21_200.txt').data(:,1));
proP = importdata('data/pro_05-21_200.txt').data(:,4);
pro = importdata('data/pro_05-21_200.txt').data(:,5);

% Prochlorococcus: 90-05
proId9 = num2str(importdata('data/proTS_90-05.txt').data(:,1));
proP9 = importdata('data/proTS_90-05.txt').data(:,4);
pro9 = importdata('data/proTS_90-05.txt').data(:,7);

%% Extract sub-ML bottle ID, pressure, bottle concentrations

% extractSMLC = 'extract sub-ML concentration'
[idSubmlChl,pSubmlChl,submlChl] = extractSMLC(chlId,chlP,chl,pMaxMld);
[idSubmlDiv,pSubmlDiv,submlDiv] = extractSMLC(divId,divP,div,pMaxMld);
[idSubmlPro,pSubmlPro,submlPro] = extractSMLC(proId,proP,pro,pMaxMld);
% [idSubmlPro9,pSubmlPro9,submlPro9] = extractSMLC(proId9,proP9,pro9,pMaxMld);

%% KS Test p-values, skewness, kurtosis: calculate

[prChl,ksChl,obsChl,skChl,kuChl,~,~,~] = ksOfLagrangian(idSubmlChl,pSubmlChl,dcm,submlChl,81);
[prDiv,ksDiv,obsDiv,skDiv,kuDiv,~,~,~] = ksOfLagrangian(idSubmlDiv,pSubmlDiv,dcm,submlDiv,69);
[prPro,ksPro,obsPro,skPro,kuPro,~,~,~] = ksOfLagrangian(idSubmlPro,pSubmlPro,dcm,submlPro,45);
% [prPro9,ksPro9,obsPro9,skPro9,kuPro9,~,~,~] = ksOfLagrangian(idSubmlPro9,pSubmlPro9,dcm,submlPro9,45);

%% KS Test p-values, skewness, kurtosis: visualise

ax1 = figure;
plotKs2(prChl,ksChl,obsChl,skChl,kuChl,prChl(1),prChl(end));
sgtitle('[chl a]');
exportgraphics(ax1,'figures/L2/ks_chla.png'); clear ax1;

ax2 = figure;
plotKs2(prDiv,ksDiv,obsDiv,skDiv,kuDiv,prDiv(1),prDiv(end));
sgtitle('[Divinyl Chl a]');
exportgraphics(ax2,'figures/L2/ks_div.png'); clear ax2;

ax3 = figure;
plotKs2(prPro,ksPro,obsPro,skPro,kuPro,prPro(1),prPro(end));
sgtitle('Prochlorococcus: 05-21');
exportgraphics(ax3,'figures/L2/ks_pro.png'); clear ax3;

% ax4 = figure;
% plotKs2(prPro9,ksPro9,obsPro9,skPro9,kuPro9,prPro9(1),prPro9(end));
% sgtitle('Prochlorococcus: 90-05');
% exportgraphics(ax4,'figures/L2/ks_pro9.png'); clear ax4;