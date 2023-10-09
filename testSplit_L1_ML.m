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
pChlIn = importdata('data/L1/hplcChla_88-21_150.txt').data(:,4);
chlIn = importdata('data/L1/hplcChla_88-21_150.txt').data(:,5);
idChlIn = importdata('data/L1/hplcChla_88-21_150.txt').data(:,1);

% 2. Extract data in ML
[idChlOut,pChlOut,chlOut] = extractMldVals(idChlIn,pChlIn,chlIn,maxMld);

% 3. Bin data
[~,pChlOutB,chlOutB,~,~] = cleanAndBin(pChlOut,chlOut,idChlOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksChl,obsChl,pChlKs,chlSk,chlKu] = ksOfBinnedCon(chlOutB,pChlOutB,10);

% 5. Plot results
ax1 = figure;
plotKs(pChlKs,ksChl,obsChl,chlSk,chlKu,0.5,20.5,true);
sgtitle('[Chl a] 88-21: L1');
exportgraphics(ax1,'figures/L1/ks_chla150.png'); clear ax1;

%% Divinyl Chlorophyll a: 88-21

% 1. Load data
pDivIn = importdata('data/L1/chlaDivi_88-21_150.txt').data(:,4);
divIn = importdata('data/L1/chlaDivi_88-21_150.txt').data(:,5);
idDivIn = importdata('data/L1/chlaDivi_88-21_150.txt').data(:,1);

% 2. Extract data in ML
[idDivOut,pDivOut,divOut] = extractMldVals(idDivIn,pDivIn,divIn,maxMld);

% 3. Bin data
[~,pDivOutB10,divOutB,~,~] = cleanAndBin(pDivOut,divOut,idDivOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksDiv,obsDiv,pDivKs,divSk,divKu] = ksOfBinnedCon(divOutB,pDivOutB10,10);

% 5. Plot results
ax2 = figure;
plotKs(pDivKs,ksDiv,obsDiv,divSk,divKu,0.5,20.5,true);
sgtitle('[Divinyl Chl a] 88-21: L1');
exportgraphics(ax2,'figures/L1/ks_divi150.png'); clear ax2;

%% Prochlorococcus: 05-21

% 1. Load data
pProIn = importdata('data/L1/pro_05-21_150.txt').data(:,4);
proIn = importdata('data/L1/pro_05-21_150.txt').data(:,5);
idProIn = importdata('data/L1/pro_05-21_150.txt').data(:,1);

% 2. Extract data in ML
[idProOut,pProOut,proOut] = extractMldVals(idProIn,pProIn,proIn,maxMld);

% 3. Bin data
[~,pProOutB10,proOutB,~,~] = cleanAndBin(pProOut,proOut,idProOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksPro,obsPro,pProKs,proSk,proKu] = ksOfBinnedCon(proOutB,pProOutB10,10);

% 5. Plot results
ax3 = figure;
plotKs(pProKs,ksPro,obsPro,proSk,proKu,0.5,20.5,true);
sgtitle('Prochlorococcus 05-21: L1');
exportgraphics(ax3,'figures/L1/ks_pro150.png'); clear ax3;

%% Prochlorococcus: 90-05

% 1. Load data
pProIn9 = importdata('data/L1/pro_90-05_150.txt').data(:,4);
proIn9 = importdata('data/L1/pro_90-05_150.txt').data(:,5);
idProIn9 = importdata('data/L1/pro_90-05_150.txt').data(:,1);

% 2. Extract data in ML
[idProOut9,pProOut9,proOut9] = extractMldVals(idProIn9,pProIn9,proIn9,maxMld);

% 3. Bin data
[~,pProOutB109,proOutB9,~,~] = cleanAndBin(pProOut9,proOut9,idProOut9');

% 4. Calculate KS p-value, skewness, kurtosis
[ksPro9,obsPro9,pProKs9,proSk9,proKu9] = ksOfBinnedCon(proOutB9,pProOutB109,10);

% 5. Plot results
ax4 = figure;
plotKs(pProKs9,ksPro9,obsPro9,proSk9,proKu9,0.5,20.5,true);
sgtitle('Prochlorococcus 90-05: L1');
exportgraphics(ax4,'figures/L1/ks_pro1509.png'); clear ax4;

%% Dissolved Oxygen: 88-21

% 1. Load data
pOxyIn = importdata('data/L1/oxy_88-21_150.txt').data(:,4);
oxyIn = importdata('data/L1/oxy_88-21_150.txt').data(:,5);
idOxyIn = importdata('data/L1/oxy_88-21_150.txt').data(:,1);

% 2. Extract data in ML
[idOxyOut,pOxyOut,oxyOut] = extractMldVals(idOxyIn,pOxyIn,oxyIn,maxMld);

% 3. Bin data
[~,pOxyOutB10,oxyOutB,~,~] = cleanAndBin(pOxyOut,oxyOut,idOxyOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksOxy,obsOxy,pOxyKs,oxySk,oxyKu] = ksOfBinnedCon(oxyOutB,pOxyOutB10,10);

% 5. Plot results
ax5 = figure;
plotKs(pOxyKs,ksOxy,obsOxy,oxySk,oxyKu,0.5,20.5,true);
sgtitle('Dissolved Oxygen 88-21: L1');
exportgraphics(ax5,'figures/L1/ks_oxy150.png'); clear ax5;

%% Dissolved Inorganic Carbon (DOC): 88-21

% 1. Load data
pDicIn = importdata('data/L1/dic_88-21_150.txt').data(:,4);
dicIn = importdata('data/L1/dic_88-21_150.txt').data(:,5);
idDicIn = importdata('data/L1/dic_88-21_150.txt').data(:,1);

% 2. Extract data in ML
[idDicOut,pDicOut,dicOut] = extractMldVals(idDicIn,pDicIn,dicIn,maxMld);

% 3. Bin data
[~,pDicOutB10,dicOutB,~,~] = cleanAndBin(pDicOut,dicOut,idDicOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksDic,obsDic,pDicKs,dicSk,dicKu] = ksOfBinnedCon(dicOutB,pDicOutB10,10);

% 5. Plot results
ax6 = figure;
plotKs(pDicKs,ksDic,obsDic,dicSk,dicKu,0.5,20.5,true);
sgtitle('Dissolved Inorganic Carbon 88-21: L1');
exportgraphics(ax6,'figures/L1/ks_dic150.png'); clear ax6;

%% pH: 92-21

% 1. Load data
pPhIn = importdata('data/L1/pH_92-21_150.txt').data(:,4);
phIn = importdata('data/L1/pH_92-21_150.txt').data(:,5);
idPhIn = importdata('data/L1/pH_92-21_150.txt').data(:,1);

% 2. Extract data in ML
[idPhOut,pPhOut,phOut] = extractMldVals(idPhIn,pPhIn,phIn,maxMld);

% 3. Bin data
[~,pPhOutB10,phOutB,~,~] = cleanAndBin(pPhOut,phOut,idPhOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksPh,obsPh,pPhKs,phSk,phKu] = ksOfBinnedCon(phOutB,pPhOutB10,10);

% 5. Plot results
ax7 = figure;
plotKs(pPhKs,ksPh,obsPh,phSk,phKu,0.5,20.5,true);
sgtitle('pH 92-21: L1');
exportgraphics(ax7,'figures/L1/ks_ph150.png'); clear ax7;

%% Alkalinity: 89-21

% 1. Load data
pAlkIn = importdata('data/L1/alk_89-21_150.txt').data(:,4);
alkIn = importdata('data/L1/alk_89-21_150.txt').data(:,5);
idAlkIn = importdata('data/L1/alk_89-21_150.txt').data(:,1);

% 2. Extract data in ML
[idAlkOut,pAlkOut,alkOut] = extractMldVals(idAlkIn,pAlkIn,alkIn,maxMld);

% 3. Bin data
[~,pAlkOutB10,alkOutB,~,~] = cleanAndBin(pAlkOut,alkOut,idAlkOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksAlk,obsAlk,pAlkKs,alkSk,alkKu] = ksOfBinnedCon(alkOutB,pAlkOutB10,10);

% 5. Plot results
ax8 = figure;
plotKs(pAlkKs,ksAlk,obsAlk,alkSk,alkKu,0.5,20.5,true);
sgtitle('Alkalinity 89-21: L1');
exportgraphics(ax8,'figures/L1/ks_alk150.png'); clear ax8;

%% Phosphate: 88-21

% 1. Load data
pPhoIn = importdata('data/L1/pho_88-21_150.txt').data(:,4);
phoIn = importdata('data/L1/pho_88-21_150.txt').data(:,5);
idPhoIn = importdata('data/L1/pho_88-21_150.txt').data(:,1);

% 2. Extract data in ML
[idPhoOut,pPhoOut,phoOut] = extractMldVals(idPhoIn,pPhoIn,phoIn,maxMld);

% 3. Bin data
[~,pPhoOutB10,phoOutB,~,~] = cleanAndBin(pPhoOut,phoOut,idPhoOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksPho,obsPho,pPhoKs,phoSk,phoKu] = ksOfBinnedCon(phoOutB,pPhoOutB10,10);

% 5. Plot results
ax9 = figure;
plotKs(pPhoKs,ksPho,obsPho,phoSk,phoKu,0.5,20.5,true);
sgtitle('Phosphate 88-21: L1');
exportgraphics(ax9,'figures/L1/ks_pho150.png'); clear ax9;

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
pSilIn = importdata('data\L1\sil_88-22_150.txt').data(:,4);
silIn = importdata('data\L1\sil_88-22_150.txt').data(:,5);
idSilIn = importdata('data\L1\sil_88-22_150.txt').data(:,1);

% 2. Extract data in ML
[idSilOut,pSilOut,silOut] = extractMldVals(idSilIn,pSilIn,silIn,maxMld);

% 3. Bin data
[~,pSilOutB10,silOutB,~,~] = cleanAndBin(pSilOut,silOut,idSilOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksSil,obsSil,pSilKs,silSk,silKu] = ksOfBinnedCon(silOutB,pSilOutB10,10);

% 5. Plot results
ax12 = figure;
plotKs(pSilKs,ksSil,obsSil,silSk,silKu);
sgtitle('Silicate 88-22: L1');
exportgraphics(ax12,'figures/L1/ks_sil150.png'); clear ax12;

%% Dissolved Organic Phosphorus: 88-01

% 1. Load data
pDopIn = importdata("data/L1/dop_88-01_150.txt").data(:,4);
dopIn = importdata("data\L1\dop_88-01_150.txt").data(:,5);
idDopIn = importdata("data\L1\dop_88-01_150.txt").data(:,1);

% 2. Extract data in ML
[idDopOut,pDopOut,dopOut] = extractMldVals(idDopIn,pDopIn,dopIn,maxMld);

% 3. Bin data
[~,pDopOutB,dopOutB,~,~] = cleanAndBin(pDopOut,dopOut,idDopOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksDop,obsDop,pDopKs,dopSk,dopKu] = ksOfBinnedCon(dopOutB,pDopOutB,10);

% 5. Plot results
ax13 = figure;
plotKs(pDopKs,ksDop,obsDop,dopSk,dopKu,0.5,20.5,true);
sgtitle('Dissolved Organic Phosphorus 88-01: L1');
exportgraphics(ax13,'figures/L1/ks_dop150.png'); clear ax13;

%% Dissolved Organic Nitrogen: 88-17

% 1. Load data
pDonIn = importdata("data\L1\don_88-17_150.txt").data(:,4);
donIn = importdata("data\L1\don_88-17_150.txt").data(:,5);
idDonIn = importdata("data\L1\don_88-17_150.txt").data(:,1);

% 2. Extract data in ML
[idDonOut,pDonOut,donOut] = extractMldVals(idDonIn,pDonIn,donIn,maxMld);

% 3. Bin data
[~,pDonOutB,donOutB,~,~] = cleanAndBin(pDonOut,donOut,idDonOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksDon,obsDon,pDonKs,donSk,donKu] = ksOfBinnedCon(donOutB,pDonOutB,10);

% 5. Plot results
ax14 = figure;
plotKs(pDonKs,ksDon,obsDon,donSk,donKu);
sgtitle('DON 88-17: L1');
exportgraphics(ax14,'figures/L1/ks_don150.png'); clear ax14;

%% Dissolved Organic Carbon (DOC): 93-17

% 1. Load data
tmp = importdata('data\L1\doc_93-17_150.txt');
pDocIn = tmp.data(:,4);
docIn = tmp.data(:,5);
idDocIn = tmp.data(:,1);
clear tmp;

% 2. Extract data in ML
[idDocOut,pDocOut,docOut] = extractMldVals(idDocIn,pDocIn,docIn,maxMld);

% 3. Bin data
[~,pDocOutB,docOutB,~,~] = cleanAndBin(pDocOut,docOut,idDocOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksDoc,obsDoc,pDocKs,docSk,docKu] = ksOfBinnedCon(docOutB,pDocOutB,10);

% 5. Plot results
ax15 = figure;
plotKs(pDocKs,ksDoc,obsDoc,docSk,docKu,0.5,20.5,true);
sgtitle('DOC 93-17: L1');
exportgraphics(ax15,'figures/L1/ks_doc150.png'); clear ax15;

%% Total Dissolved Phosphorus (TDP): 88-01

% 1. Load data
tmp = importdata('data\L1\tdp_88-01_150.txt');
pTdpIn = tmp.data(:,4);
tdpIn = tmp.data(:,5);
idTdpIn = tmp.data(:,1);
clear tmp;

% 2. Extract data in ML
[idTdpOut,pTdpOut,tdpOut] = extractMldVals(idTdpIn,pTdpIn,tdpIn,maxMld);

% 3. Bin data
[~,pTdpOutB,tdpOutB,~,~] = cleanAndBin(pTdpOut,tdpOut,idTdpOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksTdp,obsTdp,pTdpKs,tdpSk,tdpKu] = ksOfBinnedCon(tdpOutB,pTdpOutB,10);

% 5. Plot results
ax16 = figure;
plotKs(pTdpKs,ksTdp,obsTdp,tdpSk,tdpKu,0.5,20.5,true);
sgtitle('Total Dissolved Phosphorus 88-01: L1');
exportgraphics(ax16,'figures/L1/ks_tdp150.png'); clear ax16;

%% Total Dissolved Nitrogen: 88-17

% 1. Load data
tmp = importdata('data\L1\tdn_88-17_150.txt');
pTdnIn = tmp.data(:,4);
tdnIn = tmp.data(:,5);
idTdnIn = tmp.data(:,1);
clear tmp;

% 2. Extract data in ML
[idTdnOut,pTdnOut,tdnOut] = extractMldVals(idTdnIn,pTdnIn,tdnIn,maxMld);

% 3. Bin data
[~,pTdnOutB,tdnOutB,~,~] = cleanAndBin(pTdnOut,tdnOut,idTdnOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksTdn,obsTdn,pTdnKs,tdnSk,tdnKu] = ksOfBinnedCon(tdnOutB,pTdnOutB,10);

% 5. Plot results
ax = figure;
plotKs(pTdnKs,ksTdn,obsTdn,tdnSk,tdnKu);
sgtitle('Total Dissolved Nitrogen 88-17: L1');
exportgraphics(ax,'figures/L1/ks_tdn150.png'); clear ax;

%% Particulate Phosphorus: 11-21

% 1. Load data
tmp = importdata('data\L1\parp_11-21_150.txt');
pParpIn = tmp.data(:,4);
parpIn = tmp.data(:,5);
idParpIn = tmp.data(:,1);
clear tmp;

% 2. Extract data in ML
[idParpOut,pParpOut,parpOut] = extractMldVals(idParpIn,pParpIn,parpIn,maxMld);

% 3. Bin data
[~,pParpOutB,parpOutB,~,~] = cleanAndBin(pParpOut,parpOut,idParpOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksParp,obsParp,pParpKs,parpSk,parpKu] = ksOfBinnedCon(parpOutB,pParpOutB,10);

% 5. Plot results
ax = figure;
plotKs(pParpKs,ksParp,obsParp,parpSk,parpKu,0.5,20.5,true);
sgtitle('Particulate Phosphorus 11-21: L1');
exportgraphics(ax,'figures/L1/ks_parp150.png'); clear ax;

%% Particulate Nitrogen: 89-21

% 1. Load data
tmp = importdata('data\L1\parn_89-21_150.txt');
pParnIn = tmp.data(:,4);
parnIn = tmp.data(:,5);
idParnIn = tmp.data(:,1);
clear tmp;

% 2. Extract data in ML
[idParnOut,pParnOut,parnOut] = extractMldVals(idParnIn,pParnIn,parnIn,maxMld);

% 3. Bin data
[~,pParnOutB,parnOutB,~,~] = cleanAndBin(pParnOut,parnOut,idParnOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksParn,obsParn,pParnKs,parnSk,parnKu] = ksOfBinnedCon(parnOutB,pParnOutB,10);

% 5. Plot results
ax = figure;
plotKs(pParnKs,ksParn,obsParn,parnSk,parnKu,0.5,20.5,true);
sgtitle('Particulate Nitrogen 89-21: L1');
exportgraphics(ax,'figures/L1/ks_parn150.png'); clear ax;

%% Particulate Carbon: 89-21

% 1. Load data
tmp = importdata('data\L1\parc_89-21_150.txt');
pParcIn = tmp.data(:,4);
parcIn = tmp.data(:,5);
idParcIn = tmp.data(:,1);
clear tmp;

% 2. Extract data in ML
[idParcOut,pParcOut,parcOut] = extractMldVals(idParcIn,pParcIn,parcIn,maxMld);

% 3. Bin data
[~,pParcOutB,parcOutB,~,~] = cleanAndBin(pParcOut,parcOut,idParcOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ksParc,obsParc,pParcKs,parcSk,parcKu] = ksOfBinnedCon(parcOutB,pParcOutB,10);

% 5. Plot results
ax = figure;
plotKs(pParcKs,ksParc,obsParc,parcSk,parcKu,0.5,20.5,true);
sgtitle('Particulate Carbon 89-21: L1');
exportgraphics(ax,'figures/L1/ks_parc150.png'); clear ax;

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
exportgraphics(ax,'figures/L1/ks_llp150.png'); clear ax;
save("output\L1\llp.mat","p","ks","obs","Sk","Ku");

%% Low-Level Nitrogen: 89-22

% 1: Load data
tmp = importdata('data\L1\lln_89-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Low-Level Nitrogen 89-22: L1');
exportgraphics(ax,'figures/L1/ks_lln150.png'); clear ax;
save("output\L1\lln.mat","p","ks","obs","Sk","Ku");

%% Fluorometric Chlorophyll a: 89-22

% This is not of sufficient precision. Included only for reference.

% 1: Load data
tmp = importdata('data\L1\chlFlu_88-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Fluorometric Chlorophyll a 89-22: L1');
exportgraphics(ax,'figures/L1/ks_chlFlu150.png'); clear ax;
save("output\L1\chlaFluo.mat","p","ks","obs","Sk","Ku");

%% Phaeopigments: 88-22

% 1: Load data
tmp = importdata('data\L1\pheo_88-22_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Phaeopigments 88-22: L1');
exportgraphics(ax,'figures/L1/ks_phae150.png'); clear ax;
save("output\L1\phae.mat","p","ks","obs","Sk","Ku");

%% HPLC Chlorophyll C3: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\chl3_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Chl c3: L1');
exportgraphics(ax,'figures/L1/ks_chl3.png'); clear ax;
save("output\L1\chl3.mat","p","ks","obs","Sk","Ku");

%% HPLC Chlorophyll C1 + C2: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\chl12_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Chl c1 + c2: L1');
exportgraphics(ax,'figures/L1/ks_chl12.png'); clear ax;
save("output\L1\chl12.mat","p","ks","obs","Sk","Ku");

%% HPLC Chlorophyll C1 + C2 + C3: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\chl123_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Chl c1 + c2 + c3: L1');
exportgraphics(ax,'figures/L1/ks_chl123.png'); clear ax;
save("output\L1\chl123.mat","p","ks","obs","Sk","Ku");

%% Peridinin: 88-21

% One significant digit. Unreliable results!

% 1: Load data
tmp = importdata('data\L1\per_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Peridinin: L1');
exportgraphics(ax,'figures/L1/ks_per.png'); clear ax;
save("output\L1\per.mat","p","ks","obs","Sk","Ku");

%% HPLC 19' Butanoyloxyfucoxanthin: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\but19_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC 19 Butanoyloxyfucoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/ks_but19.png'); clear ax;
save("output\L1\but19.mat","p","ks","obs","Sk","Ku");

%% HPLC Fucoxanthin: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\fuco_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Fucoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/ks_fuco.png'); clear ax;
save("output\L1\fuco.mat","p","ks","obs","Sk","Ku");

%% HPLC 19' Hexanoyloxyfucoxanthin: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\hex19_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC 19'' Hexanoyloxyfucoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/ks_hex19.png'); clear ax;
save("output\L1\hex19.mat","p","ks","obs","Sk","Ku");

%% HPLC Prasinoxanthin: 88-21
% One significant digits => very unreliable!
% Also too few observations!

% 1: Load data
tmp = importdata('data\L1\prasino_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Prasinoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/ks_prasino.png'); clear ax;
save("output\L1\prasino.mat","p","ks","obs","Sk","Ku");

%% HPLC Diadinoxanthin: 88-21
% One significant digits => very unreliable!

% 1: Load data
tmp = importdata('data\L1\diadino_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Diadinoxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/ks_diadino.png'); clear ax;
save("output\L1\diadino.mat","p","ks","obs","Sk","Ku");

%% HPLC Zeaxanthin: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\zeaxan_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Zeaxanthin 88-21: L1');
exportgraphics(ax,'figures/L1/ks_zeaxan.png'); clear ax;
save("output\L1\zeaxan.mat","p","ks","obs","Sk","Ku");

%% HPLC chlorophyll b: 88-21
% Three significant digits => reliable

% 1: Load data
tmp = importdata('data\L1\chlb_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC chlorophyll b 88-21: L1');
exportgraphics(ax,'figures/L1/ks_chlb.png'); clear ax;
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
exportgraphics(ax,'figures/L1/ks_acar.png'); clear ax;
save("output\L1\acar.mat","p","ks","obs","Sk","Ku");

%% HPLC beta-Carotene: 94-21
% One or two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\bcar_94-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC beta-Carotene 94-21: L1');
exportgraphics(ax,'figures/L1/ks_bcar.png'); clear ax;
save("output\L1\bcar.mat","p","ks","obs","Sk","Ku");

%% HPLC Carotenes: 88-21
% Two significant digits => unreliable!

% 1: Load data
tmp = importdata('data\L1\caroten_88-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('HPLC Carotenes 88-21: L1');
exportgraphics(ax,'figures/L1/ks_caroten.png'); clear ax;
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

%% Prochlorococcus: 05-21
% 4 significant digits => very good!

% 1: Load data
tmp = importdata('data\L1\pbact_05-21_150.txt');

% 2-5: Extract; Bin; Calculate KS p-value, skewness, and kurtosis; and plot.
[ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);

sgtitle('Prochlorococcus 05-21: L1');
exportgraphics(ax,'figures/L1/pbact.png'); clear ax;
save("output\L1\pbact.mat","p","ks","obs","Sk","Ku");

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