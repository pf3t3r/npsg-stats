clear; clc; close all; addpath("baroneRoutines\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

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
tmp = importdata('data/L2/hplcChla_88-21_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku,sd,rV,pSubml,pV] = L2_helper(tmp,pMaxMld,dcm,[3 23]);

sgtitle('[Chl a] 88-21: L2');
exportgraphics(ax,'figures/L2/bottle/chla.png'); clear ax;
save("output\L2\chla.mat","p","ks","obs","sk","ku");
% clear ax p ks obs sk ku;

%% Prochlorococcus: 05-21

% 1. Load data
tmp = importdata('data/L2/pro_05-21_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);

sgtitle('Prochlorococcus 05-21: L2');
exportgraphics(ax,'figures/L2/bottle/pbact21.png'); clear ax;
save("output\L2\pbact21.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Prochlorococcus: 90-05

% 1. Load data
% tmp = importdata('data/L2/pro_90-05_200.txt');
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% 
% sgtitle('Prochlorococcus 90-05: L2');
% exportgraphics(ax,'figures/L2/pbact05.png'); clear ax;
% save("output\L2\pbact05.mat","p","ks","obs","sk","ku");
% clear ax p ks obs sk ku;

%% Oxygen: 88-21

% 1. Load data
tmp = importdata('data/L2/oxy_88-21_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);

sgtitle('Dissolved Oxygen 88-21: L2');
exportgraphics(ax,'figures/L2/bottle/boxy.png'); clear ax;
save("output\L2\boxy.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Dissolved Inorganic Carbon: 88-21

% 1. Load data
tmp = importdata('data/L2/dic_88-21_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);

sgtitle('Dissolved Inorganic Carbon 88-21: L2');
exportgraphics(ax,'figures/L2/bottle/dic.png'); clear ax;
save("output\L2\dic.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% pH: 92-21

% 1. Load data
tmp = importdata('data/L2/pH_92-21_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);

sgtitle('pH 92-21: L2');
exportgraphics(ax,'figures/L2/bottle/pH.png'); clear ax;
save("output\L2\pH.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Alkalinity: 89-21

% 1. Load data
tmp = importdata('data/L2/alk_89-21_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);

sgtitle('Alkalinity 89-21: L2');
exportgraphics(ax,'figures/L2/bottle/alk.png'); clear ax;
save("output\L2\alk.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Phosphate: 88-21

% 1. Load data
tmp = importdata('data/L2/pho_88-21_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);

sgtitle('Phosphate 88-21: L2');
exportgraphics(ax,'figures/L2/bottle/phos.png'); clear ax;
save("output\L2\phos.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Nitrate + Nitrite: 88-21

% 1. Load data
tmp = importdata('data/L2/nit2_88-21_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);

sgtitle('Nitrate + Nitrite 88-21: L2');
exportgraphics(ax,'figures/L2/bottle/nit2.png'); clear ax;
save("output\L2\nit2.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

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
tmp = importdata('data/L2/sil_88-22_200.txt');

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);

sgtitle('Silicate 88-22: L2');
exportgraphics(ax,'figures/L2/bottle/sil.png'); clear ax;
save("output\L2\sil.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Dissolved Organic Phosphorus (DOP): 88-01

% 1. Load data
tmp = importdata("data\L2\dop_88-01_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);

sgtitle('DOP 88-01: L2');
exportgraphics(ax,'figures/L2/bottle/dop.png'); clear ax;
save("output\L2\dop.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;


%% DON: 88-17

% Load data
tmp = importdata("data\L2\don_88-17_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('DON 88-17: L2');
exportgraphics(ax,'figures/L2/bottle/don.png'); clear ax;
save("output\L2\don.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% DOC: 93-17

% Load data
tmp = importdata("data\L2\doc_93-17_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('DOC 88-17: L2');
exportgraphics(ax,'figures/L2/bottle/doc.png'); clear ax;
save("output\L2\doc.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% TDP: 88-01

% Load data
tmp = importdata("data\L2\tdp_88-01_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('TDP 88-01: L2');
exportgraphics(ax,'figures/L2/bottle/tdp.png'); clear ax;
save("output\L2\tdp.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Total Dissolved Nitrogen (TDN): 88-17

% Load data
tmp = importdata("data\L2\tdn_88-17_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('TDN 88-17: L2');
exportgraphics(ax,'figures/L2/bottle/tdn.png'); clear ax;
save("output\L2\tdn.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Particulate Phosphorus: 11-21

% Load data
tmp = importdata("data\L2\parp_11-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('Particulate Phosphorus 11-21: L2');
exportgraphics(ax,'figures/L2/pp.png'); clear ax;
save("output\L2\pp.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Particulate Silica: 96-21

% Load data
tmp = importdata("data\L2\pars_96-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('Particulate Silica 96-21: L2');
exportgraphics(ax,'figures/L2/bottle/ps.png'); clear ax;
save("output\L2\ps.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Particulate Nitrogen: 89-21

% Load data
tmp = importdata("data\L2\parn_89-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[4 24]);
sgtitle('Particulate Nitrogen 89-21: L2');
exportgraphics(ax,'figures/L2/bottle/pn.png'); clear ax;
save("output\L2\pn.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Particulate Carbon: 89-21

% Load data
tmp = importdata("data\L2\parc_89-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[4 24]);
sgtitle('Particulate Carbon 89-21: L2');
exportgraphics(ax,'figures/L2/bottle/pc.png'); clear ax;
save("output\L2\pc.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Delta 15 N of PN

% Left as a placeholder because with the current methods I cannot deal with
% negative numbers...

%% Low-level Phosphorus: 88-22

% 1. Load data
tmp = importdata("data\L2\llp_88-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[11 31]);
sgtitle('Low-level Phosphorus 88-22: L2');
exportgraphics(ax,'figures/L2/bottle/llp.png'); clear ax;
save("output\L2\llp.mat","p","ks","obs","sk","ku");
clear ax p ks obs sk ku;

%% Low-level Nitrogen: 89-22

% 1. Load data
tmp = importdata("data\L2\lln_89-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[11 31]);
sgtitle('Low-level Nitrogen 89-22: L2');
exportgraphics(ax,'figures/L2/bottle/lln.png'); clear ax;
save("output\L2\lln.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Fluorometric Chlorophyll: 88-22

% 1. Load data
tmp = importdata("data\L2\chlFlu_88-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[11 31]);
sgtitle('Fluorometric Chlorophyll 88-22: L2');
exportgraphics(ax,'figures/L2/bottle/chlaFluo.png'); clear ax;
save("output\L2\chlaFluo.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Phaeopigments: 88-22

% 1. Load data
tmp = importdata("data\L2\pheo_88-22_200.txt");
% tmp.data = tmp.data(1:9808,:);

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[11 31]);
sgtitle('Phaeopigments 88-22: L2');
exportgraphics(ax,'figures/L2/bottle/phaeo.png'); clear ax;
save("output\L2\phaeo.mat","p","ks","obs","sk","ku");
% clear ax p ks obs sk ku;

%% HPLC Chlorophyll C3: 88-21

% 1. Load data
tmp = importdata("data\L2\chl3_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[1 21]);
sgtitle('Chlorophyll C3 (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/chl3.png'); clear ax;
save("output\L2\chl3.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Chlorophyll C1 + C2: 88-21

% 1. Load data
tmp = importdata("data\L2\chl12_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[1 21]);
sgtitle('Chlorophyll C1 + C2 (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/chl12.png'); clear ax;
save("output\L2\chl12.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Chlorophyll C1 + C2 + C3: 88-21

% 1. Load data
tmp = importdata("data\L2\chl123_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('Chlorophyll C1 + C2 + C3 (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/chl123.png'); clear ax;
save("output\L2\chl123.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Peridinin: 88-21

% ERROR: does not plot.

% % 1. Load data
% tmp = importdata("data\L2\per_88-21_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle('Peridinin (88-21): L2');
% exportgraphics(ax,'figures/L2/ks_per.png'); clear ax;
% save("output\L2\per.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% But19: 88-21

% 1. Load data
tmp = importdata("data\L2\but19_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('But19 (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/but19.png'); clear ax;
save("output\L2\but19.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Fucoxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\fuco_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Fucoxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/fuco.png'); clear ax;
save("output\L2\fuco.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC 19' Hexanoyloxyfucoxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\hex19_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC 19'' Hexanoyloxyfucoxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/hex19.png'); clear ax;
save("output\L2\hex19.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Prasinoxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\prasino_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Prasinoxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/prasino.png'); clear ax;
save("output\L2\prasino.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Diadinoxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\diadino_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Diadinoxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/diadino.png'); clear ax;
save("output\L2\diadino.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Zeaxanthin: 88-21

% 1. Load data
tmp = importdata("data\L2\zeaxan_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Zeaxanthin (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/zeaxan.png'); clear ax;
save("output\L2\zeaxan.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Chlorophyll b: 88-21

% 1. Load data
tmp = importdata("data\L2\chlb_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Chlorophyll b (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/chlb.png'); clear ax;
save("output\L2\chlb.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% 	HPLC Chlorophyll c4: 88-21
% ERROR COMPILING.

% % 1. Load data
% tmp = importdata("data\L2\chlc4_88-21_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle('HPLC Chlorophyll c4 (88-21): L2');
% exportgraphics(ax,'figures/L2/ks_chlc4.png'); clear ax;
% save("output\L2\chlc4.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC alpha-Carotene: 94-21

% 1. Load data
tmp = importdata("data\L2\acar_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC alpha-Carotene (94-21): L2');
exportgraphics(ax,'figures/L2/bottle/acar.png'); clear ax;
save("output\L2\acar.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC beta-Carotene: 94-21

% 1. Load data
tmp = importdata("data\L2\bcar_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC beta-Carotene (94-21): L2');
exportgraphics(ax,'figures/L2/bottle/bcar.png'); clear ax;
save("output\L2\bcar.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Carotenes: 88-21

% 1. Load data
tmp = importdata("data\L2\caroten_88-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Carotenes (88-21): L2');
exportgraphics(ax,'figures/L2/bottle/caroten.png'); clear ax;
save("output\L2\caroten.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC chlorophyllide a: 94-21

% 1. Load data
tmp = importdata("data\L2\chlda_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC chlorophyllide a (94-21): L2');
exportgraphics(ax,'figures/L2/bottle/chlda.png'); clear ax;
save("output\L2\chlda.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Violaxanthin: 94-21

% 1. Load data
tmp = importdata("data\L2\viol_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Violaxanthin (94-21): L2');
exportgraphics(ax,'figures/L2/bottle/viol.png'); clear ax;
save("output\L2\viol.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Lutein: 94-21

% % 1. Load data
% tmp = importdata("data\L2\lutein_94-21_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle('HPLC Lutein (94-21): L2');
% exportgraphics(ax,'figures/L2/lutein.png'); clear ax;
% save("output\L2\lutein.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Monovinyl chlorophyll a: 94-21

% 1. Load data
tmp = importdata("data\L2\mvchla_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Monovinyl chlorophyll a (94-21): L2');
exportgraphics(ax,'figures/L2/bottle/mvchla.png'); clear ax;
save("output\L2\mvchla.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% HPLC Divinyl chlorophyll a: 94-21

% 1. Load data
tmp = importdata("data\L2\dvchla_94-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('HPLC Divinyl chlorophyll a (94-21): L2');
exportgraphics(ax,'figures/L2/bottle/dvchla.png'); clear ax;
save("output\L2\dvchla.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Phycoerythrin 0.4u fraction: 00-08

% % 1. Load data
% tmp = importdata("data\L2\pe4_00-08_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle('Phycoerythrin 0.4u fraction (00-08): L2');
% exportgraphics(ax,'figures/L2/pe4.png'); clear ax;
% save("output\L2\pe4.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Phycoerythrin 5u fraction: 00-08

% % 1. Load data
% tmp = importdata("data\L2\pe5_00-08_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle('Phycoerythrin 5u fraction (00-08): L2');
% exportgraphics(ax,'figures/L2/pe5.png'); clear ax;
% save("output\L2\pe5.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Phycoerythrin 10u fraction: 00-08

% % 1. Load data
% tmp = importdata("data\L2\pe10_00-08_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle('Phycoerythrin 10u fraction (00-08): L2');
% exportgraphics(ax,'figures/L2/pe10.png'); clear ax;
% save("output\L2\pe10.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Heterotrophic Bacteria: 05-21

% 1. Load data
tmp = importdata("data\L2\hbact_05-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('Heterotrophic Bacteria (05-21): L2');
exportgraphics(ax,'figures/L2/bottle/hbact.png'); clear ax;
save("output\L2\hbact.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

% %% Prochlorococcus: 05-21
% 
% % 1. Load data
% tmp = importdata("data\L2\pbact_05-21_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle('Prochlorococcus (05-21): L2');
% exportgraphics(ax,'figures/L2/pbact.png'); clear ax;
% save("output\L2\pbact.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Synechococcus: 05-21

% % 1. Load data
% tmp = importdata("data\L2\sbact_05-21_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle('Synechococcus (05-21): L2');
% exportgraphics(ax,'figures/L2/sbact.png'); clear ax;
% save("output\L2\sbact.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Picoeukaryotes: 05-21

% 1. Load data
tmp = importdata("data\L2\ebact_05-21_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('Picoeukaryotes (05-21): L2');
exportgraphics(ax,'figures/L2/bottle/ebact.png'); clear ax;
save("output\L2\ebact.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% ATP: 88-22

% 1. Load data
tmp = importdata("data\L2\atp_88-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle('ATP (88-22): L2');
exportgraphics(ax,'figures/L2/bottle/atp.png'); clear ax;
save("output\L2\atp.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% Nitrous Oxide: 93-01

% % 1. Load data
% tmp = importdata("data\L2\n2o_93-01_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[4 14]);
% sgtitle('Nitrous Oxide (93-01): L2');
% exportgraphics(ax,'figures/L2/n2o.png'); clear ax;
% save("output\L2\n2o.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% PProd Light-12: 89-22

% 1. Load data
tmp = importdata("data\L2\l12_89-22_200.txt");

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle('PProd Light-12 (89-22): L2');
exportgraphics(ax,'figures/L2/bottle/l12.png'); clear ax;
save("output\L2\l12.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;

%% PProd Dark-12: 89-00

% % 1. Load data
% tmp = importdata("data\L2\d12_89-00_200.txt");
% 
% % Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% % and kurtosis; and plot.
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 12]);
% sgtitle('PProd Dark-12 (89-00): L2');
% exportgraphics(ax,'figures/L2/d12.png'); clear ax;
% save("output\L2\d12.mat","p","ks","obs","sk","ku"); clear ax p ks obs sk ku;
