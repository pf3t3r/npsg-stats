% Script to output L2 bottle results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 28 15]);

% Load MLD and DCM
pMaxMld = load("mldVals.mat").maxMld;
dcm = load("output/dcm.mat").dcm;  % pDcm and sigmaDcm (all casts, crn 1 - 329)

% tmpT = affix for labelling figure output.

% TEMPLATE (XX-YY)
% 1. Load data
% 2. Extract data beneath ML, centre around DCM
% 3. Calculate K-S (or A-D) p-value, Vuong LLR, skewness + kurtosis
% 4. Plot results

%% K-S
tmpT = "";

% Chlorophyll a (88-21)
tmp = importdata("data/L2/hplcChla_88-21_200.txt");
[ax,p,ks,obs,sk,ku,~,~,~,ad] = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ks");
sgtitle("[Chl a] 88-21: L2");
exportgraphics(ax,"figures/L2/bottle/chla" + tmpT + ".png");
save("output\L2\chla.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Monovinyl chlorophyll a (94-21)
tmp = importdata("data\L2\mvchla_94-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2");
exportgraphics(ax,"figures/L2/bottle/mvchla" + tmpT + ".png");
save("output\L2\mvchla.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Divinyl chlorophyll a (94-21)
tmp = importdata("data\L2\dvchla_94-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("HPLC Divinyl chlorophyll a (94-21): L2");
exportgraphics(ax,"figures/L2/bottle/dvchla" + tmpT + ".png");
save("output\L2\dvchla.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Chlorophyll b (88-21)
tmp = importdata("data\L2\chlb_88-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("HPLC Chlorophyll b (88-21): L2");
exportgraphics(ax,"figures/L2/bottle/chlb" + tmpT + ".png");
save("output\L2\chlb.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Chlorophyll C1 + C2 + C3 (88-21)
tmp = importdata("data\L2\chl123_88-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2");
exportgraphics(ax,"figures/L2/bottle/chl123" + tmpT + ".png");
save("output\L2\chl123.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC alpha-Carotene (94-21)
tmp = importdata("data\L2\acar_94-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("HPLC alpha-Carotene (94-21): L2");
exportgraphics(ax,"figures/L2/bottle/acar" + tmpT + ".png");
save("output\L2\acar.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% But19 (88-21)
tmp = importdata("data\L2\but19_88-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("But19 (88-21): L2");
exportgraphics(ax,"figures/L2/bottle/but19" + tmpT + ".png");
save("output\L2\but19.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC 19 Hexanoyloxyfucoxanthin (88-21)
tmp = importdata("data\L2\hex19_88-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2");
exportgraphics(ax,"figures/L2/bottle/hex19" + tmpT + ".png");
save("output\L2\hex19.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Zeaxanthin (88-21)
tmp = importdata("data\L2\zeaxan_88-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("HPLC Zeaxanthin (88-21): L2");
exportgraphics(ax,"figures/L2/bottle/zeaxan" + tmpT + ".png");
save("output\L2\zeaxan.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Particulate Carbon (89-21)
tmp = importdata("data\L2\parc_89-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[4 24]);
sgtitle("Particulate Carbon 89-21: L2");
exportgraphics(ax,"figures/L2/bottle/pc" + tmpT + ".png");
save("output\L2\pc.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Particulate Nitrogen (89-21)
tmp = importdata("data\L2\parn_89-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[4 24]);
sgtitle("Particulate Nitrogen 89-21: L2");
exportgraphics(ax,"figures/L2/bottle/pn" + tmpT + ".png");
save("output\L2\pn.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Nitrate + Nitrite (88-21)
tmp = importdata("data/L2/nit2_88-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("Nitrate + Nitrite 88-21: L2");
exportgraphics(ax,"figures/L2/bottle/nit2" + tmpT + ".png");
save("output\L2\nit2.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Phosphate (88-21)
tmp = importdata("data/L2/pho_88-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
sgtitle("Phosphate 88-21: L2");
exportgraphics(ax,"figures/L2/bottle/phos" + tmpT + ".png");
save("output\L2\phos.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Oxygen (88-21)
tmp = importdata("data/L2/oxy_88-21_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);
sgtitle("Dissolved Oxygen 88-21: L2");
exportgraphics(ax,"figures/L2/bottle/boxy" + tmpT + ".png");
save("output\L2\boxy.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% PProd Light-12 (89-22)
tmp = importdata("data\L2\l12_89-22_200.txt");
[ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
sgtitle("PProd Light-12 (89-22): L2");
exportgraphics(ax,"figures/L2/bottle/l12" + tmpT + ".png");
save("output\L2\l12.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

%% A-D
tmpT = "-ad";

% Chlorophyll a (88-21)
tmp = importdata("data/L2/hplcChla_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("[Chl a] 88-21: L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/chla" + tmpT + ".png");
% save("output\L2\chla.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Monovinyl chlorophyll a (94-21)
tmp = importdata("data\L2\mvchla_94-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/mvchla" + tmpT + ".png");
% save("output\L2\mvchla.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Divinyl chlorophyll a (94-21)
tmp = importdata("data\L2\dvchla_94-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/dvchla" + tmpT + ".png");
% save("output\L2\dvchla.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Chlorophyll b (88-21)
tmp = importdata("data\L2\chlb_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/chlb" + tmpT + ".png");
% save("output\L2\chlb.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Chlorophyll C1 + C2 + C3 (88-21)
tmp = importdata("data\L2\chl123_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/chl123" + tmpT + ".png");
% save("output\L2\chl123.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC alpha-Carotene (94-21)
tmp = importdata("data\L2\acar_94-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/acar" + tmpT + ".png");
% save("output\L2\acar.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% But19 (88-21)
tmp = importdata("data\L2\but19_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("But19 (88-21): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/but19" + tmpT + ".png");
% save("output\L2\but19.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC 19 Hexanoyloxyfucoxanthin (88-21)
tmp = importdata("data\L2\hex19_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/hex19" + tmpT + ".png");
% save("output\L2\hex19.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% HPLC Zeaxanthin (88-21)
tmp = importdata("data\L2\zeaxan_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("HPLC Zeaxanthin (88-21): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/zeaxan" + tmpT + ".png");
% save("output\L2\zeaxan.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Particulate Carbon (89-21)
tmp = importdata("data\L2\parc_89-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[4 24],50,4,"ad");
sgtitle("Particulate Carbon 89-21: L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/pc" + tmpT + ".png");
% save("output\L2\pc.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Particulate Nitrogen (89-21)
tmp = importdata("data\L2\parn_89-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[4 24],50,4,"ad");
sgtitle("Particulate Nitrogen 89-21: L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/pn" + tmpT + ".png");
% save("output\L2\pn.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Nitrate + Nitrite (88-21)
tmp = importdata("data/L2/nit2_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("Nitrate + Nitrite 88-21: L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/nit2" + tmpT + ".png");
% save("output\L2\nit2.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Phosphate (88-21)
tmp = importdata("data/L2/pho_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[3 23],50,4,"ad");
sgtitle("Phosphate 88-21: L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/phos" + tmpT + ".png");
% save("output\L2\phos.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% Oxygen (88-21)
tmp = importdata("data/L2/oxy_88-21_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[2 22],50,4,"ad");
sgtitle("Dissolved Oxygen 88-21: L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/boxy" + tmpT + ".png");
% save("output\L2\boxy.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

% PProd Light-12 (89-22)
tmp = importdata("data\L2\l12_89-22_200.txt");
ax = L2_helper(tmp,pMaxMld,dcm,[2 22],50,4,"ad");
sgtitle("PProd Light-12 (89-22): L2"+tmpT);
exportgraphics(ax,"figures/L2/bottle/l12" + tmpT + ".png");
% save("output\L2\l12.mat","p","ks","obs","sk","ku");
clearvars -except pMaxMld dcm tmpT;

%% K-S Unused

% % Prochlorococcus (05-21)
% tmp = importdata("data/L2/pro_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);
% sgtitle("Prochlorococcus 05-21: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/pbact21" + tmpT + ".png"); clear ax;
% save("output\L2\pbact21.mat","p","ks","obs","sk","ku");
% clear ax p ks obs sk ku;

% % Prochlorococcus (90-05)
% tmp = importdata("data/L2/pro_90-05_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("Prochlorococcus 90-05: L2");
% exportgraphics(ax,"figures/L2/pbact05.png");
% save("output\L2\pbact05.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Dissolved Inorganic Carbon (88-21)
% tmp = importdata("data/L2/dic_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);
% sgtitle("Dissolved Inorganic Carbon 88-21: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/dic" + tmpT + ".png");
% save("output\L2\dic.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % pH (92-21)
% tmp = importdata("data/L2/pH_92-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);
% sgtitle("pH 92-21: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/pH" + tmpT + ".png");
% save("output\L2\pH.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Alkalinity (89-21)
% tmp = importdata("data/L2/alk_89-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 22]);
% sgtitle("Alkalinity 89-21: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/alk" + tmpT + ".png");
% save("output\L2\alk.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Silicate (88-22)
% tmp = importdata("data/L2/sil_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("Silicate 88-22: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/sil" + tmpT + ".png");
% save("output\L2\sil.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Dissolved Organic Phosphorus (DOP) (88-01)
% tmp = importdata("data\L2\dop_88-01_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("DOP 88-01: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/dop" + tmpT + ".png");
% save("output\L2\dop.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % DON (88-17)
% tmp = importdata("data\L2\don_88-17_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("DON 88-17: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/don" + tmpT + ".png"); clear ax;
% save("output\L2\don.mat","p","ks","obs","sk","ku");
% clear ax p ks obs sk ku tmpT;

% % DOC (93-17)
% tmp = importdata("data\L2\doc_93-17_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("DOC 88-17: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/doc" + tmpT + ".png");
% save("output\L2\doc.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % TDP (88-01)
% tmp = importdata("data\L2\tdp_88-01_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("TDP 88-01: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/tdp" + tmpT + ".png");
% save("output\L2\tdp.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Total Dissolved Nitrogen (TDN) (88-17)
% tmp = importdata("data\L2\tdn_88-17_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("TDN 88-17: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/tdn" + tmpT + ".png");
% save("output\L2\tdn.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Particulate Phosphorus (11-21)
% tmp = importdata("data\L2\parp_11-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("Particulate Phosphorus 11-21: L2");
% exportgraphics(ax,"figures/L2/pp" + tmpT + ".png");
% save("output\L2\pp.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Particulate Silica (96-21)
% tmp = importdata("data\L2\pars_96-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("Particulate Silica 96-21: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/ps" + tmpT + ".png");
% save("output\L2\ps.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Low-level Phosphorus (88-22)
% tmp = importdata("data\L2\llp_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[11 31]);
% sgtitle("Low-level Phosphorus 88-22: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/llp" + tmpT + ".png");
% save("output\L2\llp.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Low-level Nitrogen (89-22)
% tmp = importdata("data\L2\lln_89-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[11 31]);
% sgtitle("Low-level Nitrogen 89-22: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/lln" + tmpT + ".png");
% save("output\L2\lln.mat","p","ks","obs","sk","ku"); clearvars -except pMaxMld dcm tmpT;

% % Fluorometric Chlorophyll (88-22)
% tmp = importdata("data\L2\chlFlu_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[11 31]);
% sgtitle("Fluorometric Chlorophyll 88-22: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/chlaFluo" + tmpT + ".png");
% save("output\L2\chlaFluo.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Phaeopigments (88-22)
% tmp = importdata("data\L2\pheo_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[11 31]);
% sgtitle("Phaeopigments 88-22: L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/phaeo" + tmpT + ".png");
% save("output\L2\phaeo.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC Chlorophyll C3 (88-21)
% tmp = importdata("data\L2\chl3_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[1 21]);
% sgtitle("Chlorophyll C3 (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/chl3" + tmpT + ".png");
% save("output\L2\chl3.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC Chlorophyll C1 + C2 (88-21)
% tmp = importdata("data\L2\chl12_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[1 21]);
% sgtitle("Chlorophyll C1 + C2 (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/chl12" + tmpT + ".png");
% save("output\L2\chl12.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC Fucoxanthin (88-21)
% tmp = importdata("data\L2\fuco_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("HPLC Fucoxanthin (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/fuco" + tmpT + ".png");
% save("output\L2\fuco.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC Prasinoxanthin (88-21)
% tmp = importdata("data\L2\prasino_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("HPLC Prasinoxanthin (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/prasino" + tmpT + ".png");
% save("output\L2\prasino.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC Diadinoxanthin (88-21)
% tmp = importdata("data\L2\diadino_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("HPLC Diadinoxanthin (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/diadino" + tmpT + ".png");
% save("output\L2\diadino.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC beta-Carotene (94-21)
% tmp = importdata("data\L2\bcar_94-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("HPLC beta-Carotene (94-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/bcar" + tmpT + ".png");
% save("output\L2\bcar.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC Carotenes (88-21)
% tmp = importdata("data\L2\caroten_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("HPLC Carotenes (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/caroten" + tmpT + ".png");
% save("output\L2\caroten.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC chlorophyllide a (94-21)
% tmp = importdata("data\L2\chlda_94-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("HPLC chlorophyllide a (94-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/chlda" + tmpT + ".png");
% save("output\L2\chlda.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % HPLC Violaxanthin (94-21)
% tmp = importdata("data\L2\viol_94-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("HPLC Violaxanthin (94-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/viol" + tmpT + ".png");
% save("output\L2\viol.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% HPLC Lutein (94-21)
% tmp = importdata("data\L2\lutein_94-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("HPLC Lutein (94-21): L2");
% exportgraphics(ax,"figures/L2/lutein.png");
% save("output\L2\lutein.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Phycoerythrin 0.4u fraction (00-08)
% tmp = importdata("data\L2\pe4_00-08_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("Phycoerythrin 0.4u fraction (00-08): L2");
% exportgraphics(ax,"figures/L2/pe4.png");
% save("output\L2\pe4.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Phycoerythrin 5u fraction (00-08)
% tmp = importdata("data\L2\pe5_00-08_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("Phycoerythrin 5u fraction (00-08): L2");
% exportgraphics(ax,"figures/L2/pe5.png");
% save("output\L2\pe5.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Phycoerythrin 10u fraction (00-08)
% tmp = importdata("data\L2\pe10_00-08_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("Phycoerythrin 10u fraction (00-08): L2");
% exportgraphics(ax,"figures/L2/pe10.png");
% save("output\L2\pe10.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Heterotrophic Bacteria (05-21)
% tmp = importdata("data\L2\hbact_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("Heterotrophic Bacteria (05-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/hbact" + tmpT + ".png");
% save("output\L2\hbact.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Prochlorococcus (05-21)
% tmp = importdata("data\L2\pbact_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("Prochlorococcus (05-21): L2");
% exportgraphics(ax,"figures/L2/pbact.png");
% save("output\L2\pbact.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Synechococcus (05-21)
% tmp = importdata("data\L2\sbact_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("Synechococcus (05-21): L2");
% exportgraphics(ax,"figures/L2/sbact.png");
% save("output\L2\sbact.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Picoeukaryotes (05-21)
% tmp = importdata("data\L2\ebact_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm);
% sgtitle("Picoeukaryotes (05-21): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/ebact" + tmpT + ".png");
% save("output\L2\ebact.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % ATP (88-22)
% tmp = importdata("data\L2\atp_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[3 23]);
% sgtitle("ATP (88-22): L2");
% exportgraphics(ax,"figures/L2/bottle/notUsed/atp" + tmpT + ".png");
% save("output\L2\atp.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % Nitrous Oxide (93-01)
% tmp = importdata("data\L2\n2o_93-01_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[4 14]);
% sgtitle("Nitrous Oxide (93-01): L2");
% exportgraphics(ax,"figures/L2/n2o.png");
% save("output\L2\n2o.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;

% % PProd Dark-12 (89-00)
% tmp = importdata("data\L2\d12_89-00_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,pMaxMld,dcm,[2 12]);
% sgtitle("PProd Dark-12 (89-00): L2");
% exportgraphics(ax,"figures/L2/d12.png");
% save("output\L2\d12.mat","p","ks","obs","sk","ku");
% clearvars -except pMaxMld dcm tmpT;