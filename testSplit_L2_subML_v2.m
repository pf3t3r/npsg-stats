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

%% Particulate Silica: 96-21

% Load data
tmp = importdata("data\L2\pars_96-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,44);
sgtitle('Particulate Silica 96-21: L2');
exportgraphics(ax,'figures/L2/ks_pars.png'); clear ax;
save("output\L2\pars.mat","p","ks","obs","sk","ku");
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

%% HPLC Fucoxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\fuco_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC Fucoxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/ks_fuco.png'); clear ax;
save("output\L2\fuco.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC 19' Hexanoyloxyfucoxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\hex19_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC 19'' Hexanoyloxyfucoxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/ks_hex19.png'); clear ax;
save("output\L2\hex19.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Prasinoxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\prasino_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,41);
sgtitle('HPLC Prasinoxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/ks_prasino.png'); clear ax;
save("output\L2\prasino.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Diadinoxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\diadino_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC Diadinoxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/ks_diadino.png'); clear ax;
save("output\L2\diadino.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Zeaxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\zeaxan_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC Zeaxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/ks_zeaxan.png'); clear ax;
save("output\L2\zeaxan.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Chlorophyll b: 88-21

% 1. Load data
tmp = importdata("data\L2\chlb_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC Chlorophyll b (88-21): L2');
exportgraphics(ax,'figures/L2/ks_chlb.png'); clear ax;
save("output\L2\chlb.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% 	HPLC Chlorophyll c4: 88-21
% ERROR COMPILING.

% 1. Load data
tmp = importdata("data\L2\chlc4_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC Chlorophyll c4 (88-21): L2');
exportgraphics(ax,'figures/L2/ks_chlc4.png'); clear ax;
save("output\L2\chlc4.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC alpha-Carotene: 94-21

% 1. Load data
tmp = importdata("data\L2\acar_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC alpha-Carotene (94-21): L2');
exportgraphics(ax,'figures/L2/ks_acar.png'); clear ax;
save("output\L2\acar.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC beta-Carotene: 94-21

% 1. Load data
tmp = importdata("data\L2\bcar_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC beta-Carotene (94-21): L2');
exportgraphics(ax,'figures/L2/ks_bcar.png'); clear ax;
save("output\L2\bcar.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Carotenes: 88-21

% 1. Load data
tmp = importdata("data\L2\caroten_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC Carotenes (88-21): L2');
exportgraphics(ax,'figures/L2/ks_caroten.png'); clear ax;
save("output\L2\caroten.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC chlorophyllide a: 94-21

% 1. Load data
tmp = importdata("data\L2\chlda_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,35);
sgtitle('HPLC chlorophyllide a (94-21): L2');
exportgraphics(ax,'figures/L2/ks_chlda.png'); clear ax;
save("output\L2\chlda.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Violaxanthin: 94-21

% 1. Load data
tmp = importdata("data\L2\viol_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,40);
sgtitle('HPLC Violaxanthin (94-21): L2');
exportgraphics(ax,'figures/L2/viol.png'); clear ax;
save("output\L2\viol.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Lutein: 94-21

% 1. Load data
tmp = importdata("data\L2\lutein_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,29);
sgtitle('HPLC Lutein (94-21): L2');
exportgraphics(ax,'figures/L2/lutein.png'); clear ax;
save("output\L2\lutein.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Monovinyl chlorophyll a: 94-21

% 1. Load data
tmp = importdata("data\L2\mvchla_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC Monovinyl chlorophyll a (94-21): L2');
exportgraphics(ax,'figures/L2/mvchla.png'); clear ax;
save("output\L2\mvchla.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Divinyl chlorophyll a: 94-21

% 1. Load data
tmp = importdata("data\L2\dvchla_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('HPLC Divinyl chlorophyll a (94-21): L2');
exportgraphics(ax,'figures/L2/dvchla.png'); clear ax;
save("output\L2\dvchla.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Phycoerythrin 0.4u fraction: 00-08

% 1. Load data
tmp = importdata("data\L2\pe4_00-08_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('Phycoerythrin 0.4u fraction (00-08): L2');
exportgraphics(ax,'figures/L2/pe4.png'); clear ax;
save("output\L2\pe4.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Phycoerythrin 5u fraction: 00-08

% 1. Load data
tmp = importdata("data\L2\pe5_00-08_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('Phycoerythrin 5u fraction (00-08): L2');
exportgraphics(ax,'figures/L2/pe5.png'); clear ax;
save("output\L2\pe5.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Phycoerythrin 10u fraction: 00-08

% 1. Load data
tmp = importdata("data\L2\pe10_00-08_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('Phycoerythrin 10u fraction (00-08): L2');
exportgraphics(ax,'figures/L2/pe10.png'); clear ax;
save("output\L2\pe10.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Heterotrophic Bacteria: 05-21

% 1. Load data
tmp = importdata("data\L2\hbact_05-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,42);
sgtitle('Heterotrophic Bacteria (05-21): L2');
exportgraphics(ax,'figures/L2/hbact.png'); clear ax;
save("output\L2\hbact.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Prochlorococcus: 05-21

% 1. Load data
tmp = importdata("data\L2\pbact_05-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,42);
sgtitle('Prochlorococcus (05-21): L2');
exportgraphics(ax,'figures/L2/pbact.png'); clear ax;
save("output\L2\pbact.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Synechococcus: 05-21

% 1. Load data
tmp = importdata("data\L2\sbact_05-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,29);
sgtitle('Synechococcus (05-21): L2');
exportgraphics(ax,'figures/L2/sbact.png'); clear ax;
save("output\L2\sbact.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Picoeukaryotes: 05-21

% 1. Load data
tmp = importdata("data\L2\ebact_05-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,42);
sgtitle('Picoeukaryotes (05-21): L2');
exportgraphics(ax,'figures/L2/ebact.png'); clear ax;
save("output\L2\ebact.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% ATP: 88-22

% 1. Load data
tmp = importdata("data\L2\atp_88-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,69);
sgtitle('ATP (88-22): L2');
exportgraphics(ax,'figures/L2/atp.png'); clear ax;
save("output\L2\atp.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Nitrous Oxide: 93-01

% 1. Load data
tmp = importdata("data\L2\n2o_93-01_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('Nitrous Oxide (93-01): L2');
exportgraphics(ax,'figures/L2/n2o.png'); clear ax;
save("output\L2\n2o.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% PProd Light-12: 89-22

% 1. Load data
tmp = importdata("data\L2\l12_89-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,42);
sgtitle('PProd Light-12 (89-22): L2');
exportgraphics(ax,'figures/L2/l12.png'); clear ax;
save("output\L2\l12.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% PProd Dark-12: 89-00

% 1. Load data
tmp = importdata("data\L2\d12_89-00_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,42);
sgtitle('PProd Dark-12 (89-00): L2');
exportgraphics(ax,'figures/L2/d12.png'); clear ax;
save("output\L2\d12.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;