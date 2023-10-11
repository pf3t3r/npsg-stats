clear; clc; close all; addpath("baroneRoutines\");

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

%% Extract Maximum Mixed Layer Depth (per cruise) 'maxMld'

ctdData = importdata('datafiles\ctd_iso_ALL.mat').ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData i;

%% README
% How is this script organised?
% At the highest level we have each variable in 'blocks'.
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

%% Chlorophyll a: 88-21

% 1. Load data
tmp = importdata('data/L1/hplcChla_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('[Chl a] 88-21: L1');
exportgraphics(ax,'figures/L1/chla.png'); clear ax;
save("output\L1\chla.mat","p","ks","obs","Sk","Ku");

%% Prochlorococcus: 05-21

% 1. Load data
tmp = importdata('data/L1/pro_05-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Prochlorococcus 05-21: L1');
exportgraphics(ax,'figures/L1/pbact21.png'); clear ax;
save("output\L1\pbact21.mat","p","ks","obs","Sk","Ku");

%% Prochlorococcus: 90-05

% 1. Load data
tmp = importdata('data/L1/pro_90-05_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Prochlorococcus 90-05: L1');
exportgraphics(ax,'figures/L1/pbact05.png'); clear ax;
save("output\L1\pbact05.mat","p","ks","obs","Sk","Ku");

%% Bottle Dissolved Oxygen: 88-21

% 1. Load data
tmp = importdata('data/L1/oxy_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Dissolved Oxygen 88-21: L1');
exportgraphics(ax,'figures/L1/boxy.png'); clear ax;
save("output\L1\boxy.mat","p","ks","obs","Sk","Ku");

%% Dissolved Inorganic Carbon (DOC): 88-21

% 1. Load data
tmp = importdata('data/L1/dic_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Dissolved Inorganic Carbon 88-21: L1');
exportgraphics(ax,'figures/L1/dic.png'); clear ax;
save("output\L1\dic.mat","p","ks","obs","Sk","Ku");

%% pH: 92-21

% 1. Load data
tmp = importdata('data/L1/pH_92-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('pH 92-21: L1');
exportgraphics(ax,'figures/L1/pH.png'); clear ax;
save("output\L1\pH.mat","p","ks","obs","Sk","Ku");

%% Alkalinity: 89-21

% 1. Load data
tmp = importdata('data/L1/alk_89-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Alkalinity 89-21: L1');
exportgraphics(ax,'figures/L1/alk.png'); clear ax;
save("output\L1\alk.mat","p","ks","obs","Sk","Ku");

%% Phosphate: 88-21

% 1. Load data
tmp = importdata('data/L1/pho_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Phosphate 88-21: L1');
exportgraphics(ax,'figures/L1/phos.png'); clear ax;
save("output\L1\phos.mat","p","ks","obs","Sk","Ku");

%% Nitrate + Nitrite: 88-21

% 1. Load data
pNit2In = importdata('data/L1/nit2_88-21_150.txt').data(:,4);
nit2In = importdata('data/L1/nit2_88-21_150.txt').data(:,5);
idNit2In = importdata('data/L1/nit2_88-21_150.txt').data(:,1);

% 2. Extract data in ML
[idNit2Out,pNit2Out,nit2Out] = extractMldVals(idNit2In,pNit2In,nit2In,maxMld);

% 3. Bin data
[~,pNit2OutB,nit2OutB,~,~] = cleanAndBin(pNit2Out,nit2Out,idNit2Out');

% 4. Calculate KS p-value, skewness, kurtosis
[ksNit2,obsNit2,pNit2Ks,nit2Sk,nit2Ku] = ksOfBinnedCon(nit2OutB,pNit2OutB,10);

% 5. Plot results
% % Broken: no pressures output so can't plot.
% ax10 = figure; % Nitrate + Nitrite: 88-21
% plotKs(pNit2Ks,ksNit2,obsNit2,nit2Sk,nit2Ku,0.5,20.5,true);
% sgtitle('Nitrate + Nitrite 88-21: Mixed Layer');
% exportgraphics(ax10,'figures/L1/ks_nit2_150.png'); clear ax10;

%% Nitrite: 89-94

% 1. Load data
pNitIn = importdata('data/L1/nit_89-95_150.txt').data(:,4);
nitIn = importdata('data/L1/nit_89-95_150.txt').data(:,5);
idNitIn = importdata('data/L1/nit_89-95_150.txt').data(:,1);

% 2. Extract data in ML
% Nitrite: 89-94 ''' breaks
% [idNitOut,pNitOut,nitOut] = extractMldVals(idNitIn,pNitIn,nitIn,maxMld);

% 3. Bin data
% Nitrite: 89-94
% [~,pNitOutB10,nitOutB,~,~] = cleanAndBin(pNitOut,nitOut,idNitOut');

% 4. Calculate KS p-value, skewness, kurtosis
% Nitrite: 89-94
% [ksNit,obsNit,pNitKs,nitSk,nitKu] = ksOfBinnedCon(nitOutB,pNitOutB10,10);

% 5. Plot results
% ax11 = figure; % Nitrite: 89-94
% plotKs(pNitKs,ksNit,obsNit,nitSk,nitKu,0.5,20.5,true);
% sgtitle('Nitrite 89-94: Mixed Layer');
% exportgraphics(ax11,'figures/L1/ks_nit150.png'); clear ax11;

%% Silicate: 88-22

% 1. Load data
tmp = importdata('data\L1\sil_88-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Silicate 88-22: L1');
exportgraphics(ax,'figures/L1/sil.png'); clear ax;
save("output\L1\sil.mat","p","ks","obs","Sk","Ku");

%% Dissolved Organic Phosphorus: 88-01

% 1. Load data
tmp = importdata("data/L1/dop_88-01_150.txt");

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Dissolved Organic Phosphorus 88-01: L1');
exportgraphics(ax,'figures/L1/dop.png'); clear ax;
save("output\L1\dop.mat","p","ks","obs","Sk","Ku");

%% Dissolved Organic Nitrogen: 88-17

% 1. Load data
tmp = importdata("data\L1\don_88-17_150.txt");

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('DON 88-17: L1');
exportgraphics(ax,'figures/L1/don.png'); clear ax;
save("output\L1\don.mat","p","ks","obs","Sk","Ku");

%% Dissolved Organic Carbon (DOC): 93-17

% 1. Load data
tmp = importdata('data\L1\doc_93-17_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('DOC 93-17: L1');
exportgraphics(ax,'figures/L1/doc.png'); clear ax;
save("output\L1\doc.mat","p","ks","obs","Sk","Ku");

%% Total Dissolved Phosphorus (TDP): 88-01

% 1. Load data
tmp = importdata('data\L1\tdp_88-01_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Total Dissolved Phosphorus 88-01: L1');
exportgraphics(ax,'figures/L1/tdp.png'); clear ax;
save("output\L1\tdp.mat","p","ks","obs","Sk","Ku");

%% Total Dissolved Nitrogen: 88-17

% 1. Load data
tmp = importdata('data\L1\tdn_88-17_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Total Dissolved Nitrogen 88-17: L1');
exportgraphics(ax,'figures/L1/tdn.png'); clear ax;
save("output\L1\tdn.mat","p","ks","obs","Sk","Ku");

%% Particulate Phosphorus: 11-21

% 1. Load data
tmp = importdata('data\L1\parp_11-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Particulate Phosphorus 11-21: L1');
exportgraphics(ax,'figures/L1/pp.png'); clear ax;
save("output\L1\pp.mat","p","ks","obs","Sk","Ku");

%% Particulate Nitrogen: 89-21

% 1. Load data
tmp = importdata('data\L1\parn_89-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Particulate Nitrogen 89-21: L1');
exportgraphics(ax,'figures/L1/pn.png'); clear ax;
save("output\L1\pn.mat","p","ks","obs","Sk","Ku");

%% Particulate Carbon: 89-21

% 1. Load data
tmp = importdata('data\L1\parc_89-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Particulate Carbon 89-21: L1');
exportgraphics(ax,'figures/L1/pc.png'); clear ax;
save("output\L1\pc.mat","p","ks","obs","Sk","Ku");

%% delta15N of PN: 00-04

% This DOES NOT work because of the negative values.

% % 1. Load data
% tmp = importdata('data\L1\p15n_00-04_150.txt');
% pP15nIn = tmp.data(:,4);
% p15nIn = tmp.data(:,5);
% idP15nIn = tmp.data(:,1);
% clear tmp;
% 
% % 2. Extract data in ML
% [idP15nOut,pP15nOut,p15nOut] = extractMldVals(idP15nIn,pP15nIn,p15nIn,maxMld);
% 
% % 3. Bin data
% [~,pP15nOutB,p15nOutB,~,~] = cleanAndBin(pP15nOut,p15nOut,idP15nOut');
% 
% % 4. Calculate KS p-value, skewness, kurtosis
% [ksP15n,obsP15n,pP15n,p15nSk,p15nKu] = ksOfBinnedCon(p15nOutB,pP15nOutB,10,39);
% 
% % 5. Plot results
% ax = figure;
% plotKs(pP15n,ksP15n,obsP15n,p15nSk,p15nKu,0.5,20.5,true,39);
% sgtitle('delta 15 N of PN 00-04: L1');
% exportgraphics(ax,'figures/L1/ks_p15n150.png'); clear ax;

%% Low-Level Phosphorus: 88-22

% 1: Load data
tmp = importdata('data\L1\llp_88-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Low-Level Phosphorus 88-22: L1');
exportgraphics(ax,'figures/L1/llp.png'); clear ax;
save("output\L1\llp.mat","p","ks","obs","Sk","Ku");

%% Low-Level Nitrogen: 89-22

% 1: Load data
tmp = importdata('data\L1\lln_89-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Low-Level Nitrogen 89-22: L1');
exportgraphics(ax,'figures/L1/lln.png'); clear ax;
save("output\L1\lln.mat","p","ks","obs","Sk","Ku");

%% Fluorometric Chlorophyll a: 89-22

% This is not of sufficient precision. Included only for reference.

% 1: Load data
tmp = importdata('data\L1\chlFlu_88-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Fluorometric Chlorophyll a 89-22: L1');
exportgraphics(ax,'figures/L1/chlaFluo.png'); clear ax;
save("output\L1\chlaFluo.mat","p","ks","obs","Sk","Ku");

%% Phaeopigments: 88-22

% 1: Load data
tmp = importdata('data\L1\pheo_88-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Phaeopigments 88-22: L1');
exportgraphics(ax,'figures/L1/phaeo.png'); clear ax;
save("output\L1\phaeo.mat","p","ks","obs","Sk","Ku");

%% HPLC Chlorophyll C3: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\chl3_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Chl c3: L1');
exportgraphics(ax,'figures/L1/chl3.png'); clear ax;
save("output\L1\chl3.mat","p","ks","obs","Sk","Ku");

%% HPLC Chlorophyll C1 + C2: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\chl12_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Chl c1 + c2: L1');
exportgraphics(ax,'figures/L1/chl12.png'); clear ax;
save("output\L1\chl12.mat","p","ks","obs","Sk","Ku");

%% HPLC Chlorophyll C1 + C2 + C3: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\chl123_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Chl c1 + c2 + c3: L1');
exportgraphics(ax,'figures/L1/chl123.png'); clear ax;
save("output\L1\chl123.mat","p","ks","obs","Sk","Ku");

%% Peridinin: 88-21

% One significant digit. Unreliable results!

% 1: Load data
tmp = importdata('data\L1\per_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Peridinin: L1');
exportgraphics(ax,'figures/L1/perid.png'); clear ax;
save("output\L1\perid.mat","p","ks","obs","Sk","Ku");

%% HPLC 19' Butanoyloxyfucoxanthin: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\but19_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC 19 Butanoyloxyfucoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/but19.png'); clear ax;
save("output\L1\but19.mat","p","ks","obs","Sk","Ku");

%% HPLC Fucoxanthin: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\fuco_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Fucoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/fuco.png'); clear ax;
save("output\L1\fuco.mat","p","ks","obs","Sk","Ku");

%% HPLC 19' Hexanoyloxyfucoxanthin: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\hex19_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC 19'' Hexanoyloxyfucoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/hex19.png'); clear ax;
save("output\L1\hex19.mat","p","ks","obs","Sk","Ku");

%% HPLC Prasinoxanthin: 88-21
% One significant digits => very unreliable!
% Also too few observations!

% 1: Load data
tmp = importdata('data\L1\prasino_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Prasinoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/prasino.png'); clear ax;
save("output\L1\prasino.mat","p","ks","obs","Sk","Ku");

%% HPLC Diadinoxanthin: 88-21
% One significant digits => very unreliable!

% 1: Load data
tmp = importdata('data\L1\diadino_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Diadinoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/diadino.png'); clear ax;
save("output\L1\diadino.mat","p","ks","obs","Sk","Ku");

%% HPLC Zeaxanthin: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\zeaxan_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Zeaxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/zeaxan.png'); clear ax;
save("output\L1\zeaxan.mat","p","ks","obs","Sk","Ku");

%% HPLC chlorophyll b: 88-21
% Three significant digits => reliable

% 1: Load data
tmp = importdata('data\L1\chlb_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC chlorophyll b 88-21: L1');
exportgraphics(ax,'figures/L1/chlb.png'); clear ax;
save("output\L1\chlb.mat","p","ks","obs","Sk","Ku");

%% HPLC chlorophyll c4: 88-21
% ERROR COMPILING
% Two significant digits => unreliable!

% % 1: Load data
% tmp = importdata('data\L1\chlc4_88-21_150.txt');
% 
% % 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% 
% sgtitle('HPLC chlorophyll c4 88-21: L1');
% exportgraphics(ax,'figures/L1/ks_chlc4.png'); clear ax;
% save("output\L1\chlc4.mat","p","ks","obs","Sk","Ku");

%% HPLC alpha-Carotene: 94-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\acar_94-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC alpha-Carotene 94-21: L1');
exportgraphics(ax,'figures/L1/acar.png'); clear ax;
save("output\L1\acar.mat","p","ks","obs","Sk","Ku");

%% HPLC beta-Carotene: 94-21
% One or two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\bcar_94-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC beta-Carotene 94-21: L1');
exportgraphics(ax,'figures/L1/bcar.png'); clear ax;
save("output\L1\bcar.mat","p","ks","obs","Sk","Ku");

%% HPLC Carotenes: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\caroten_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Carotenes 88-21: L1');
exportgraphics(ax,'figures/L1/caroten.png'); clear ax;
save("output\L1\caroten.mat","p","ks","obs","Sk","Ku");

%% HPLC Chlorophyllide a: 94-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\chlda_94-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Chlorophyllide a 94-21: L1');
exportgraphics(ax,'figures/L1/chlda.png'); clear ax;
save("output\L1\chlda.mat","p","ks","obs","Sk","Ku");

%% HPLC Violaxanthin: 94-21
% One significant digits => very unreliable!
% ERROR: parameters must be positive.

% % 1: Load data
% tmp = importdata('data\L1\viol_94-21_150.txt');
% 
% % 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% 
% sgtitle('HPLC Violaxanthin 94-21: L1');
% exportgraphics(ax,'figures/L1/viol.png'); clear ax;
% save("output\L1\viol.mat","p","ks","obs","Sk","Ku");

%% HPLC Lutein: 94-21
% One significant digits => very unreliable!

% 1: Load data
tmp = importdata('data\L1\lutein_94-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Lutein 94-21: L1');
exportgraphics(ax,'figures/L1/lutein.png'); clear ax;
save("output\L1\lutein.mat","p","ks","obs","Sk","Ku");

%% HPLC Monovinyl chlorophyll a: 94-21
% Three significant digits => good!

% 1: Load data
tmp = importdata('data\L1\mvchla_94-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Monovinyl chlorophyll a 94-21: L1');
exportgraphics(ax,'figures/L1/mvchla.png'); clear ax;
save("output\L1\mvchla.mat","p","ks","obs","Sk","Ku");

%% HPLC Divinyl chlorophyll a: 94-21
% Three significant digits => good!

% 1: Load data
tmp = importdata('data\L1\dvchla_94-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Divinyl chlorophyll a 94-21: L1');
exportgraphics(ax,'figures/L1/dvchla.png'); clear ax;
save("output\L1\dvchla.mat","p","ks","obs","Sk","Ku");

%% Phycoerythrin 0.4u fraction: 00-08
% Four significant digits => very good!

% 1: Load data
tmp = importdata('data\L1\pe4_00-08_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Phycoerythrin 0.4u fraction 00-08: L1');
exportgraphics(ax,'figures/L1/pe4.png'); clear ax;
save("output\L1\pe4.mat","p","ks","obs","Sk","Ku");

%% Phycoerythrin 5u fraction: 00-08
% Three significant digits => good!

% 1: Load data
tmp = importdata('data\L1\pe5_00-08_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Phycoerythrin 5u fraction 00-08: L1');
exportgraphics(ax,'figures/L1/pe5.png'); clear ax;
save("output\L1\pe5.mat","p","ks","obs","Sk","Ku");

%% Phycoerythrin 10u fraction: 00-08
% 4 significant digits => very good!

% 1: Load data
tmp = importdata('data\L1\pe10_00-08_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Phycoerythrin 10u fraction 00-08: L1');
exportgraphics(ax,'figures/L1/pe10.png'); clear ax;
save("output\L1\pe10.mat","p","ks","obs","Sk","Ku");

%% Heterotrophic Bacteria: 05-21
% 4 significant digits => very good!

% 1: Load data
tmp = importdata('data\L1\hbact_05-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Heterotrophic Bacteria 05-21: L1');
exportgraphics(ax,'figures/L1/hbact.png'); clear ax;
save("output\L1\hbact.mat","p","ks","obs","Sk","Ku");

% %% Prochlorococcus: 05-21
% % 4 significant digits => very good!
% 
% % 1: Load data
% tmp = importdata('data\L1\pbact_05-21_150.txt');
% 
% % 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% 
% sgtitle('Prochlorococcus 05-21: L1');
% exportgraphics(ax,'figures/L1/pbact.png'); clear ax;
% save("output\L1\pbact.mat","p","ks","obs","Sk","Ku");

%% Synechococcus: 05-21
% Two significant digits => may be unreliable!

% 1: Load data
tmp = importdata('data\L1\sbact_05-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Synechococcus 05-21: L1');
exportgraphics(ax,'figures/L1/sbact.png'); clear ax;
save("output\L1\sbact.mat","p","ks","obs","Sk","Ku");

%% PicoEukaryotes: 05-21
% 2 significant digits => may be unreliable!

% 1: Load data
tmp = importdata('data\L1\ebact_05-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Picoeukaryotes 05-21: L1');
exportgraphics(ax,'figures/L1/ebact.png'); clear ax;
save("output\L1\ebact.mat","p","ks","obs","Sk","Ku");

%% ATP: 88-22
% Four significant digits => very good!

% 1: Load data
tmp = importdata('data\L1\atp_88-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('ATP 88-22: L1');
exportgraphics(ax,'figures/L1/atp.png'); clear ax;
save("output\L1\atp.mat","p","ks","obs","Sk","Ku");

%% Nitrous Oxide: 93-01
% 4 significant digits => very good!

% 1: Load data
tmp = importdata('data\L1\n2o_93-01_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Nitrous Oxide 93-01: L1');
exportgraphics(ax,'figures/L1/n2o.png'); clear ax;
save("output\L1\n2o.mat","p","ks","obs","Sk","Ku");

%% PProd Light-12: 89-22
% 4 significant digits => very good!

% 1: Load data
tmp = importdata('data\L1\l12_89-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('PProd Light-12 89-22: L1');
exportgraphics(ax,'figures/L1/l12.png'); clear ax;
save("output\L1\l12.mat","p","ks","obs","Sk","Ku");

%% PProd Dark-12: 89-00
% Three significant digits => good!

% 1: Load data
tmp = importdata('data\L1\d12_89-00_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('PProd Dark-12 89-00: L1');
exportgraphics(ax,'figures/L1/d12.png'); clear ax;
save("output\L1\d12.mat","p","ks","obs","Sk","Ku");

%% Visualise ML Extraction

% pressure = pChlOut;
% concentration = chlOut;
% xLabel = 'Chl a [ng/l]';
% figTitle = 'Chl a: 88-21';
% 
% figure; % Chlorophyll a
% scatter(concentration,pressure);
% grid on;
% set(gca,'YDir','reverse'); ylabel('Pressure [dbar]');
% xlabel(xLabel);
% title(figTitle);

%% Save new data

save mldVals.mat maxMld;