clear; clc; close all; addpath("baroneRoutines\");

%% Load MLD and DCM
pMaxMld = load('mldVals.mat').maxMld;
dcm = load("dcm.mat").dcm;  % pDcm and sigmaDcm (all casts, crn 1 - 329)

%% TEMPLATE: XX-YY

% 1. Load data

% 2. Extract data beneath ML, centre around DCM

% 3. Calculate KS p-value, skewness, kurtosis

% 4. Plot results

%% Chlorophyll a: 88-21

% 1. Load data
chlId = num2str(importdata('data/L2/hplcChla_88-21_200.txt').data(:,1));
chlP = importdata('data/L2/hplcChla_88-21_200.txt').data(:,4);
chl = importdata('data/L2/hplcChla_88-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlChl,pSubmlChl,submlChl] = extractSMLC(chlId,chlP,chl,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prChl,ksChl,obsChl,skChl,kuChl,~,~,~] = ksOfLagrangian(idSubmlChl,pSubmlChl,dcm,submlChl,81);

% 4. Plot results
ax1 = figure;
plotKs2(prChl,ksChl,obsChl,skChl,kuChl,prChl(1),prChl(end),81);
sgtitle('[Chl a] 88-21: L2');
exportgraphics(ax1,'figures/L2/ks_chla.png'); clear ax1;

%% Divinyl Chlorophyll a: 88-21

% 1. Load data
divId = num2str(importdata('data/L2/diviChla_88-21_200.txt').data(:,1));
divP = importdata('data/L2/diviChla_88-21_200.txt').data(:,4);
div = importdata('data/L2/diviChla_88-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlDiv,pSubmlDiv,submlDiv] = extractSMLC(divId,divP,div,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prDiv,ksDiv,obsDiv,skDiv,kuDiv,~,~,~] = ksOfLagrangian(idSubmlDiv,pSubmlDiv,dcm,submlDiv,69);

% 4. Plot results
ax2 = figure;
plotKs2(prDiv,ksDiv,obsDiv,skDiv,kuDiv,prDiv(1),prDiv(end),69);
sgtitle('[Divinyl Chl a] 88-21: L2');
exportgraphics(ax2,'figures/L2/ks_div.png'); clear ax2;

%% Prochlorococcus: 05-21

% 1. Load data
proId = num2str(importdata('data/L2/pro_05-21_200.txt').data(:,1));
proP = importdata('data/L2/pro_05-21_200.txt').data(:,4);
pro = importdata('data/L2/pro_05-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlPro,pSubmlPro,submlPro] = extractSMLC(proId,proP,pro,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prPro,ksPro,obsPro,skPro,kuPro,~,~,~] = ksOfLagrangian(idSubmlPro,pSubmlPro,dcm,submlPro,42);

% 4. Plot results
ax3 = figure;
plotKs2(prPro,ksPro,obsPro,skPro,kuPro,prPro(1),prPro(end),42);
sgtitle('Prochlorococcus 05-21: L2');
exportgraphics(ax3,'figures/L2/ks_pro.png'); clear ax3;

%% Prochlorococcus: 90-05

% 1. Load data
proId9 = num2str(importdata('data/L2/pro_90-05_200.txt').data(:,1));
proP9 = importdata('data/L2/pro_90-05_200.txt').data(:,4);
pro9 = importdata('data/L2/pro_90-05_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlPro9,pSubmlPro9,submlPro9] = extractSMLC(proId9,proP9,pro9,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prPro9,ksPro9,obsPro9,skPro9,kuPro9,~,~,~] = ksOfLagrangian(idSubmlPro9,pSubmlPro9,dcm,submlPro9,30);

% 4. Plot results
ax4 = figure;
plotKs2(prPro9,ksPro9,obsPro9,skPro9,kuPro9,prPro9(1),prPro9(end),30);
sgtitle('Prochlorococcus 90-05: L2');
exportgraphics(ax4,'figures/L2/ks_pro9.png'); clear ax4;

%% Oxygen: 88-21

% 1. Load data
oxyId = num2str(importdata('data/L2/oxy_88-21_200.txt').data(:,1));
oxyP = importdata('data/L2/oxy_88-21_200.txt').data(:,4);
oxy = importdata('data/L2/oxy_88-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlOxy,pSubmlOxy,submlOxy] = extractSMLC(oxyId,oxyP,oxy,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prOxy,ksOxy,obsOxy,skOxy,kuOxy,~,~,~] = ksOfLagrangian(idSubmlOxy,pSubmlOxy,dcm,submlOxy,122);

% 4. Plot results
ax5 = figure;
plotKs2(prOxy,ksOxy,obsOxy,skOxy,kuOxy,prOxy(1),prOxy(end),122);
sgtitle('Dissolved Oxygen 88-21: L2');
exportgraphics(ax5,'figures/L2/ks_oxy.png'); clear ax5;

%% Dissolved Inorganic Carbon: 88-21

% 1. Load data
dicId = num2str(importdata('data/L2/dic_88-21_200.txt').data(:,1));
dicP = importdata('data/L2/dic_88-21_200.txt').data(:,4);
dic = importdata('data/L2/dic_88-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlDic,pSubmlDic,submlDic] = extractSMLC(dicId,dicP,dic,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prDic,ksDic,obsDic,skDic,kuDic,~,~,~] = ksOfLagrangian(idSubmlDic,pSubmlDic,dcm,submlDic,89);

% 4. Plot results
ax6 = figure;
plotKs2(prDic,ksDic,obsDic,skDic,kuDic,prDic(1),prDic(end),89);
sgtitle('Dissolved Inorganic Carbon 88-21: L2');
exportgraphics(ax6,'figures/L2/ks_dic.png'); clear ax6;

%% pH: 92-21

% 1. Load data
phId = num2str(importdata('data/L2/pH_92-21_200.txt').data(:,1));
phP = importdata('data/L2/pH_92-21_200.txt').data(:,4);
pH = importdata('data/L2/pH_92-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlPh,pSubmlPh,submlPh] = extractSMLC(phId,phP,pH,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prPh,ksPh,obsPh,skPh,kuPh,~,~,~] = ksOfLagrangian(idSubmlPh,pSubmlPh,dcm,submlPh,59);

% 4. Plot results
ax7 = figure;
plotKs2(prPh,ksPh,obsPh,skPh,kuPh,prPh(1),prPh(end),59);
sgtitle('pH 92-21: L2');
exportgraphics(ax7,'figures/L2/ks_pH.png'); clear ax7;

%% Alkalinity: 89-21

% 1. Load data
alkId = num2str(importdata('data/L2/alk_89-21_200.txt').data(:,1));
alkP = importdata('data/L2/alk_89-21_200.txt').data(:,4);
alk = importdata('data/L2/alk_89-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlAlk,pSubmlAlk,submlAlk] = extractSMLC(alkId,alkP,alk,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prAlk,ksAlk,obsAlk,skAlk,kuAlk,~,~,~] = ksOfLagrangian(idSubmlAlk,pSubmlAlk,dcm,submlAlk,83);

% 4. Plot results
ax8 = figure;
plotKs2(prAlk,ksAlk,obsAlk,skAlk,kuAlk,prAlk(1),prAlk(end),83);
sgtitle('Alkalinity 89-21: L2');
exportgraphics(ax8,'figures/L2/ks_alk.png'); clear ax8;

%% Phosphate: 88-21

% 1. Load data
phoId = num2str(importdata('data/L2/pho_88-21_200.txt').data(:,1));
phoP = importdata('data/L2/pho_88-21_200.txt').data(:,4);
pho = importdata('data/L2/pho_88-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlPho,pSubmlPho,submlPho] = extractSMLC(phoId,phoP,pho,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prPho,ksPho,obsPho,skPho,kuPho,~,~,~] = ksOfLagrangian(idSubmlPho,pSubmlPho,dcm,submlPho);

% 4. Plot results
ax9 = figure;
plotKs2(prPho,ksPho,obsPho,skPho,kuPho,prPho(1),prPho(end));
sgtitle('Phosphate 88-21: L2');
exportgraphics(ax9,'figures/L2/ks_pho.png'); clear ax9;

%% Nitrate + Nitrite: 88-21

% 1. Load data
nit2Id = num2str(importdata('data/L2/nit2_88-21_200.txt').data(:,1));
nit2P = importdata('data/L2/nit2_88-21_200.txt').data(:,4);
nit2 = importdata('data/L2/nit2_88-21_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlNit2,pSubmlNit2,submlNit2] = extractSMLC(nit2Id,nit2P,nit2,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prNit2,ksNit2,obsNit2,skNit2,kuNit2,~,~,~] = ksOfLagrangian(idSubmlNit2,pSubmlNit2,dcm,submlNit2);

% 4. Plot results
ax10 = figure;
plotKs2(prNit2,ksNit2,obsNit2,skNit2,kuNit2,prNit2(1),prNit2(end));
sgtitle('Nitrate + Nitrite 88-21: L2');
exportgraphics(ax10,'figures/L2/ks_nit2.png'); clear ax10;

%% Nitrite: 89-94

% 1. Load data
nitId = num2str(importdata('data/L2/nit_89-94_200.txt').data(:,1));
nitP = importdata('data/L2/nit_89-94_200.txt').data(:,4);
nit = importdata('data/L2/nit_89-94_200.txt').data(:,5);

% % 2. Extract data beneath ML, centre around DCM
% [idSubmlNit,pSubmlNit,submlNit] = extractSMLC(nitId,nitP,nit,pMaxMld);

% % 3. Calculate KS p-value, skewness, kurtosis
% [prNit,ksNit,obsNit,skNit,kuNit,~,~,~] = ksOfLagrangian(idSubmlNit,pSubmlNit,dcm,submlNit);

% % 4. Plot results
% ax11 = figure;
% plotKs2(prNit,ksNit,obsNit,skNit,kuNit,prNit(1),prNit(end));
% sgtitle('Nitrite 89-94: L2');
% exportgraphics(ax11,'figures/L2/ks_nit.png'); clear ax11;

%% Silicate: 88-22

% 1. Load data
silId = num2str(importdata('data/L2/sil_88-22_200.txt').data(:,1));
silP = importdata('data/L2/sil_88-22_200.txt').data(:,4);
sil = importdata('data/L2/sil_88-22_200.txt').data(:,5);

% 2. Extract data beneath ML, centre around DCM
[idSubmlSil,pSubmlSil,submlSil] = extractSMLC(silId,silP,sil,pMaxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[prSil,ksSil,obsSil,skSil,kuSil,~,~,~] = ksOfLagrangian(idSubmlSil,pSubmlSil,dcm,submlSil);

% 4. Plot results
ax12 = figure; % Silicate: 88-22
plotKs2(prSil,ksSil,obsSil,skSil,kuSil,prSil(1),prSil(end));
sgtitle('Silicate 88-22: L2');
exportgraphics(ax12,'figures/L2/ks_sil.png'); clear ax12;

%% Dissolved Organic Phosphorus (DOP): 88-01

% 1. Load data
tmp = importdata("data\L2\dop_88-01_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,44);
sgtitle('DOP 88-01: L2');
exportgraphics(ax,'figures/L2/ks_dop.png'); clear ax;
save("output\L2\dop.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;


%% DON: 88-17

% Load data
tmp = importdata("data\L2\don_88-17_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,44);
sgtitle('DON 88-17: L2');
exportgraphics(ax,'figures/L2/ks_don.png'); clear ax;
save("output\L2\don.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% DOC: 93-17

% Load data
tmp = importdata("data\L2\doc_93-17_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,65);
sgtitle('DOC 88-17: L2');
exportgraphics(ax,'figures/L2/ks_doc.png'); clear ax;
save("output\L2\doc.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% TDP: 88-01

% Load data
tmp = importdata("data\L2\tdp_88-01_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,44);
sgtitle('TDP 88-01: L2');
exportgraphics(ax,'figures/L2/ks_tdp.png'); clear ax;
save("output\L2\tdp.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Total Dissolved Nitrogen (TDN): 88-17

% Load data
tmp = importdata("data\L2\tdn_88-17_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,56);
sgtitle('TDN 88-17: L2');
exportgraphics(ax,'figures/L2/ks_tdn.png'); clear ax;
save("output\L2\tdn.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Particulate Phosphorus: 11-21

% Load data
tmp = importdata("data\L2\parp_11-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,20);
sgtitle('Particulate Phosphorus 11-21: L2');
exportgraphics(ax,'figures/L2/ks_parp.png'); clear ax;
save("output\L2\parp.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Particulate Nitrogen: 89-21

% Load data
tmp = importdata("data\L2\parn_89-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,55);
sgtitle('Particulate Nitrogen 89-21: L2');
exportgraphics(ax,'figures/L2/ks_parn.png'); clear ax;
save("output\L2\parn.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Particulate Carbon: 89-21

% Load data
tmp = importdata("data\L2\parc_89-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,55);
sgtitle('Particulate Carbon 89-21: L2');
exportgraphics(ax,'figures/L2/ks_parc.png'); clear ax;
save("output\L2\parc.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Delta 15 N of PN

% Left as a placeholder because with the current methods I cannot deal with
% negative numbers...

%% Low-level Phosphorus: 88-22

% 1. Load data
tmp = importdata("data\L2\llp_88-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,100);
sgtitle('Low-level Phosphorus 88-22: L2');
exportgraphics(ax,'figures/L2/ks_llp.png'); clear ax;
save("output\L2\llp.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Low-level Nitrogen: 89-22

% 1. Load data
tmp = importdata("data\L2\lln_89-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,78);
sgtitle('Low-level Nitrogen 89-22: L2');
exportgraphics(ax,'figures/L2/ks_lln.png'); clear ax;
save("output\L2\lln.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Fluorometric Chlorophyll: 88-22

% 1. Load data
tmp = importdata("data\L2\chlFlu_88-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,100);
sgtitle('Fluorometric Chlorophyll 88-22: L2');
exportgraphics(ax,'figures/L2/ks_chlaFluo.png'); clear ax;
save("output\L2\chlaFluo.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Phaeopigments: 88-22

% 1. Load data
tmp = importdata("data\L2\pheo_88-22_200.txt");
% tmp.data = tmp.data(1:9808,:);

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('Phaeopigments 88-22: L2');
exportgraphics(ax,'figures/L2/ks_phae.png'); clear ax;
save("output\L2\phae.mat","p","ks","obs","sk","ku");
% clear ax p ks obs sk ku;

%% HPLC Chlorophyll C3: 88-21

% 1. Load data
tmp = importdata("data\L2\chl3_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,21);
sgtitle('Chlorophyll C3 (88-21): L2');
exportgraphics(ax,'figures/L2/ks_chl3.png'); clear ax;
save("output\L2\chl3.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Chlorophyll C1 + C2: 88-21

% 1. Load data
tmp = importdata("data\L2\chl12_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,28);
sgtitle('Chlorophyll C1 + C2 (88-21): L2');
exportgraphics(ax,'figures/L2/ks_chl12.png'); clear ax;
save("output\L2\chl12.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Chlorophyll C1 + C2 + C3: 88-21

% 1. Load data
tmp = importdata("data\L2\chl123_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,94);
sgtitle('Chlorophyll C1 + C2 + C3 (88-21): L2');
exportgraphics(ax,'figures/L2/ks_chl123.png'); clear ax;
save("output\L2\chl123.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Peridinin: 88-21

% 1. Load data
tmp = importdata("data\L2\per_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,94);
sgtitle('Peridinin (88-21): L2');
exportgraphics(ax,'figures/L2/ks_per.png'); clear ax;
save("output\L2\per.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% But19: 88-21

% 1. Load data
tmp = importdata("data\L2\but19_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,98);
sgtitle('But19 (88-21): L2');
exportgraphics(ax,'figures/L2/ks_but19.png'); clear ax;
save("output\L2\but19.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;
