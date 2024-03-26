% Script to output L1 bottle results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\");
addpath("func\");

set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 28 15]);

% Label for figures
% tmpT = "NL";
tmpT = "_S";

%% Extract Maximum Mixed Layer Depth (per cruise) "maxMld"

ctdData = importdata("datafiles\ctd_iso_ALL.mat").ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData i;

%% Load 3rd and 4th Moments Bias
% unused

% % For now use "Case 3" (parameters based on estimation for chla at 100
% % dbar)
% tmpK = importdata("output\skku\kurtBiasC3.mat");
% tmpS = importdata("output\skku\skewBiasC3.mat");
% 
% uN = [tmpK.n1'-tmpK.sn1(:,1) tmpK.sn1(:,2)-tmpK.n1' tmpS.mn1'-tmpS.ssn1(:,1) tmpS.ssn1(:,2)-tmpS.mn1'];
% uL = [tmpK.l1'-tmpK.sl1(:,1) tmpK.sl1(:,2)-tmpK.l1' tmpS.ml1'-tmpS.ssl1(:,1) tmpS.ssl1(:,2)-tmpS.ml1'];
% uG = [tmpK.g1'-tmpK.pg1(:,1) tmpK.pg1(:,2)-tmpK.g1' tmpS.g2'-tmpS.pg2(:,1) tmpS.pg2(:,2)-tmpS.g2'];
% uW = [tmpK.w1'-tmpK.pw1(:,1) tmpK.pw1(:,2)-tmpK.w1' tmpS.w2'-tmpS.pw2(:,1) tmpS.pw2(:,2)-tmpS.w2'];
% 
% unc = [uN uL uG uW];
%% README
% How is this script organised?
% At the highest level we have each variable in "blocks".
% Within these blocks we 
% (1) load the data, (2) extract that data which lies in the mixed layer,
% (3) bin this data, (4) calculate the KS p-value, skewness, and kurtosis
% for this binned mixed-layer data, and finally (5) plot (and save?) the
% results.
% THE TWO SECTIONS BEFORE THE README MUST BE RUN IN ORDER TO LOAD THE MLD,
% THEN EACH INDIVIDUAL BLOCK CAN BE RUN.
%% TEMPLATE: XX-YY

% 1. Load data

% 2. Extract data in ML

% 3. Bin data

% 4. Calculate KS p-value, skewness, kurtosis

% 5. Plot results

%% K-S Analysis.

% chl-a (88-21)
tmp = importdata("data/L1/hplcChla_88-21_150.txt");
[ax,p,ks,obs,Sk,Ku,~,~] = L1_helper(tmp,maxMld);
sgtitle("[Chl a] 88-21: L1");
exportgraphics(ax,"figures/L1/bottle/chla" + tmpT + ".png");
save("output\L1\chla.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% HPLC Monovinyl chlorophyll a (94-21)
tmp = importdata("data\L1\mvchla_94-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1");
exportgraphics(ax,"figures/L1/bottle/mvchla" + tmpT + ".png");
save("output\L1\mvchla.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% HPLC Divinyl chlorophyll a (94-21)
tmp = importdata("data\L1\dvchla_94-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("HPLC Divinyl chlorophyll a 94-21: L1");
exportgraphics(ax,"figures/L1/bottle/dvchla" + tmpT + ".png");
save("output\L1\dvchla.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% HPLC chlorophyll b (88-21)
tmp = importdata("data\L1\chlb_88-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("HPLC chlorophyll b 88-21: L1");
exportgraphics(ax,"figures/L1/bottle/chlb" + tmpT + ".png");
save("output\L1\chlb.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% HPLC Chlorophyll C1 + C2 + C3 (88-21)
tmp = importdata("data\L1\chl123_88-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("HPLC Chl c1 + c2 + c3: L1");
exportgraphics(ax,"figures/L1/bottle/chl123" + tmpT + ".png");
save("output\L1\chl123.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% HPLC alpha-Carotene (94-21)
tmp = importdata("data\L1\acar_94-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("HPLC alpha-Carotene 94-21: L1");
exportgraphics(ax,"figures/L1/bottle/acar" + tmpT + ".png");
save("output\L1\acar.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% HPLC 19' Butanoyloxyfucoxanthin (88-21)
tmp = importdata("data\L1\but19_88-21_150.txt");
[ax,p,ks,obs,Sk,Ku,~,~,~] = L1_helper(tmp,maxMld);
sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1");
exportgraphics(ax,"figures/L1/bottle/but19" + tmpT + ".png");
save("output\L1\but19.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% HPLC 19" Hexanoyloxyfucoxanthin (88-21)
tmp = importdata("data\L1\hex19_88-21_150.txt");
[ax,p,ks,obs,Sk,Ku,~,~,~] = L1_helper(tmp,maxMld);
sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1");
exportgraphics(ax,"figures/L1/bottle/hex19" + tmpT + ".png");
save("output\L1\hex19.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% HPLC Zeaxanthin (88-21)
tmp = importdata("data\L1\zeaxan_88-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("HPLC Zeaxanthin 88-21: L1");
exportgraphics(ax,"figures/L1/bottle/zeaxan" + tmpT + ".png");
save("output\L1\zeaxan.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% Particulate Carbon (89-21)
tmp = importdata("data\L1\parc_89-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("Particulate Carbon 89-21: L1");
exportgraphics(ax,"figures/L1/bottle/pc" + tmpT + ".png");
save("output\L1\pc.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% Particulate Nitrogen (89-21)
tmp = importdata("data\L1\parn_89-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("Particulate Nitrogen 89-21: L1");
exportgraphics(ax,"figures/L1/bottle/pn" + tmpT + ".png");
save("output\L1\pn.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% Low-Level Phosphorus (88-22)
tmp = importdata("data\L1\llp_88-22_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("Low-Level Phosphorus 88-22: L1");
exportgraphics(ax,"figures/L1/bottle/llp" + tmpT + ".png");
save("output\L1\llp.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% Low-Level Nitrogen (89-22)
tmp = importdata("data\L1\lln_89-22_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("Low-Level Nitrogen 89-22: L1");
exportgraphics(ax,"figures/L1/bottle/lln" + tmpT + ".png");
save("output\L1\lln.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% Bottle Dissolved Oxygen (88-21)
tmp = importdata("data/L1/oxy_88-21_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("Dissolved Oxygen 88-21: L1");
exportgraphics(ax,"figures/L1/bottle/boxy" + tmpT + ".png");
save("output\L1\boxy.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

% PProd Light-12 (89-22) (4 significant digits => very good!)
tmp = importdata("data\L1\l12_89-22_150.txt");
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
sgtitle("PProd Light-12 89-22: L1");
exportgraphics(ax,"figures/L1/bottle/l12" + tmpT + ".png");
save("output\L1\l12.mat","p","ks","obs","Sk","Ku");
clearvars -except tmpT maxMld;

%% Visualise ML Extraction

% pressure = pChlOut;
% concentration = chlOut;
% xLabel = "Chl a [ng/l]";
% figTitle = "Chl a: 88-21";
% 
% figure; % Chlorophyll a
% scatter(concentration,pressure);
% grid on;
% set(gca,"YDir","reverse"); ylabel("Pressure [dbar]");
% xlabel(xLabel);
% title(figTitle);

%% Save new data

save mldVals.mat maxMld;

%% A-D Analysis
tmpT = "-ad";

% HPLC chl-a
tmp = importdata("data/L1/hplcChla_88-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("[Chl a] 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/chla" + tmpT + ".png");
clearvars -except tmpT maxMld;

% mvchla
tmp = importdata("data/L1/mvchla_94-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("[Monovinyl Chl a] 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/mvchla" + tmpT + ".png");
clearvars -except tmpT maxMld;

% dvchla
tmp = importdata("data/L1/dvchla_94-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("[Divinyl Chl a] 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/dvchla" + tmpT + ".png");
clearvars -except tmpT maxMld;

% chl-b
tmp = importdata("data/L1/chlb_88-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("[Chl b] 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/chlb" + tmpT + ".png");
clearvars -except tmpT maxMld;

% chl-c123
tmp = importdata("data/L1/chl123_88-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("[Chl-c 1 + 2 + 3] 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/chlc123" + tmpT + ".png");
clearvars -except tmpT maxMld;

% alpha-caro
tmp = importdata("data/L1/acar_94-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("alpha-carotene 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/chlb" + tmpT + ".png");
clearvars -except tmpT maxMld;

% but-19
tmp = importdata("data/L1/but19_88-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("But-19 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/but19" + tmpT + ".png");
clearvars -except tmpT maxMld;

% hex-19
tmp = importdata("data/L1/hex19_88-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("Hex-19 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/hex19" + tmpT + ".png");
clearvars -except tmpT maxMld;

% zeax
tmp = importdata("data/L1/zeaxan_88-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("Zeaxanthin 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/zeax" + tmpT + ".png");
clearvars -except tmpT maxMld;

% pc
tmp = importdata("data/L1/parc_89-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("Particulate Carbon 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/pc" + tmpT + ".png");
clearvars -except tmpT maxMld;

% pn
tmp = importdata("data/L1/parn_89-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("Particulate Nitrogen 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/pn" + tmpT + ".png");
clearvars -except tmpT maxMld;

% llp l1
tmp = importdata("data/L1/llp_88-22_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("Low-Level Phosphorus 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/llp" + tmpT + ".png");
clearvars -except tmpT maxMld;

% lln l1
tmp = importdata("data/L1/lln_89-22_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("Low-Level Nitrogen 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/lln" + tmpT + ".png");
clearvars -except tmpT maxMld;

% phosphate
tmp = importdata("data\L1\pho_88-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("Phosphate 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/lln" + tmpT + ".png");
clearvars -except tmpT maxMld;

% do
tmp = importdata("data/L1/oxy_88-21_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("Dissolved Oxygen 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/boxy" + tmpT + ".png");
clearvars -except tmpT maxMld;

% l-12 pp
tmp = importdata("data/L1/l12_89-22_150.txt");
ax = L1_helper(tmp,maxMld,50,4,"ad");
sgtitle("L-12 PP: L1"+tmpT);
exportgraphics(ax,"figures/L1/bottle/l12" + tmpT + ".png");
clearvars -except tmpT maxMld;

%% Parameters not used for current analysis

% % Dissolved Organic Nitrogen (88-17)
% tmp = importdata("data\L1\don_88-17_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("DON 88-17: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/don" + tmpT + ".png");
% save("output\L1\don.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Dissolved Organic Carbon (DOC): 93-17
% tmp = importdata("data\L1\doc_93-17_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("DOC 93-17: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/doc" + tmpT + ".png");
% save("output\L1\doc.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Total Dissolved Phosphorus (TDP): 88-01
% tmp = importdata("data\L1\tdp_88-01_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Total Dissolved Phosphorus 88-01: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/tdp" + tmpT + ".png");
% save("output\L1\tdp.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Total Dissolved Nitrogen: 88-17
% tmp = importdata("data\L1\tdn_88-17_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Total Dissolved Nitrogen 88-17: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/tdn" + tmpT + ".png");
% save("output\L1\tdn.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Particulate Phosphorus: 11-21
% tmp = importdata("data\L1\parp_11-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Particulate Phosphorus 11-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/pp" + tmpT + ".png");
% save("output\L1\pp.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Fluorometric Chlorophyll a (89-22)
% tmp = importdata("data\L1\chlFlu_88-22_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Fluorometric Chlorophyll a 89-22: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/chlaFluo" + tmpT + ".png");
% save("output\L1\chlaFluo.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Phaeopigments (88-22)
% tmp = importdata("data\L1\pheo_88-22_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phaeopigments 88-22: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/phaeo" + tmpT + ".png");
% save("output\L1\phaeo.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % HPLC Chlorophyll C3 (88-21)
% tmp = importdata("data\L1\chl3_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Chl c3: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/chl3" + tmpT + ".png");
% save("output\L1\chl3.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % HPLC Chlorophyll C1 + C2: 88-21
% tmp = importdata("data\L1\chl12_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld); 
% sgtitle("HPLC Chl c1 + c2: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/chl12" + tmpT + ".png");
% save("output\L1\chl12.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Peridinin (88-21)
% tmp = importdata("data\L1\per_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Peridinin: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/perid" + tmpT + ".png");
% save("output\L1\perid.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % HPLC Fucoxanthin (88-21)
% tmp = importdata("data\L1\fuco_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Fucoxanthin 88-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/fuco" + tmpT + ".png"); clear ax;
% save("output\L1\fuco.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % HPLC Diadinoxanthin (88-21)
% tmp = importdata("data\L1\diadino_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Diadinoxanthin 88-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/diadino" + tmpT + ".png");
% save("output\L1\diadino.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % HPLC beta-Carotene (94-21)
% tmp = importdata("data\L1\bcar_94-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC beta-Carotene 94-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/bcar" + tmpT + ".png");
% save("output\L1\bcar.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % HPLC Carotenes (88-21)
% tmp = importdata("data\L1\caroten_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Carotenes 88-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/caroten" + tmpT + ".png"); clear ax;
% save("output\L1\caroten.mat","p","ks","obs","Sk","Ku");

% % HPLC Chlorophyllide a (94-21)
% tmp = importdata("data\L1\chlda_94-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Chlorophyllide a 94-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/chlda" + tmpT + ".png");
% save("output\L1\chlda.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % HPLC Lutein (94-21)
% tmp = importdata("data\L1\lutein_94-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Lutein 94-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/lutein" + tmpT + ".png");
% save("output\L1\lutein.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Phycoerythrin 0.4u fraction (00-08)
% tmp = importdata("data\L1\pe4_00-08_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phycoerythrin 0.4u fraction 00-08: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/pe4" + tmpT + ".png");
% save("output\L1\pe4.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Phycoerythrin 5u fraction (00-08) (Three significant digits => good!)
% tmp = importdata("data\L1\pe5_00-08_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phycoerythrin 5u fraction 00-08: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/pe5" + tmpT + ".png");
% save("output\L1\pe5.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Phycoerythrin 10u fraction (00-08) (4 significant digits => very good!)
% tmp = importdata("data\L1\pe10_00-08_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phycoerythrin 10u fraction 00-08: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/pe10" + tmpT + ".png");
% save("output\L1\pe10.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Heterotrophic Bacteria (05-21) (4 significant digits => very good!)
% tmp = importdata("data\L1\hbact_05-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Heterotrophic Bacteria 05-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/hbact" + tmpT + ".png");
% save("output\L1\hbact.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Prochlorococcus (05-21) (4 significant digits => very good!)
% tmp = importdata("data\L1\pbact_05-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Prochlorococcus 05-21: L1");
% exportgraphics(ax,"figures/L1/pbact.png"); clear ax;
% save("output\L1\pbact.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % Synechococcus (05-21) (Two significant digits => may be unreliable!)
% tmp = importdata("data\L1\sbact_05-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Synechococcus 05-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/sbact" + tmpT + ".png");
% save("output\L1\sbact.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % PicoEukaryotes: 05-21 (2 significant digits => may be unreliable!)
% tmp = importdata("data\L1\ebact_05-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Picoeukaryotes 05-21: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/ebact" + tmpT + ".png");
% save("output\L1\ebact.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % ATP (88-22) (Four significant digits => very good!)
% tmp = importdata("data\L1\atp_88-22_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("ATP 88-22: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/atp" + tmpT + ".png");
% save("output\L1\atp.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% % PProd Dark-12 (89-00) (Three significant digits => good!)
% tmp = importdata("data\L1\d12_89-00_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("PProd Dark-12 89-00: L1");
% exportgraphics(ax,"figures/L1/bottle/notUsed/d12" + tmpT + ".png");
% save("output\L1\d12.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;

% Dissolved Organic Phosphate 
% tmp = importdata("data\L1\pho_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phosphate 88-21: L1");
% exportgraphics(ax,"figures/L1/bottle/phos" + tmpT + ".png");
% save("output\L1\phos.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld;