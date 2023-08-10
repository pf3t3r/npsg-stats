clear; clc; close all; addpath("baroneRoutines\");

%% Extract Maximum Mixed Layer Depth (per cruise) 'maxMld'

ctdData = importdata('datafiles\ctd_iso_ALL.mat').ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData;

%% Load Data

% Chlorophyll a
pChlIn = importdata('data/L1/hplcChla_88-21_150.txt').data(:,4);
chlIn = importdata('data/L1/hplcChla_88-21_150.txt').data(:,5);
idChlIn = importdata('data/L1/hplcChla_88-21_150.txt').data(:,1);

% Divinyl Chl a
pDivIn = importdata('data/L1/chlaDivi_88-21_150.txt').data(:,4);
divIn = importdata('data/L1/chlaDivi_88-21_150.txt').data(:,5);
idDivIn = importdata('data/L1/chlaDivi_88-21_150.txt').data(:,1);

% Prochlorococcus 05-21
pProIn = importdata('data/L1/pro_05-21_150.txt').data(:,4);
proIn = importdata('data/L1/pro_05-21_150.txt').data(:,5);
idProIn = importdata('data/L1/pro_05-21_150.txt').data(:,1);

% Prochlorococcus 90-05
pProIn9 = importdata('data/L1/pro_90-05_150.txt').data(:,4);
proIn9 = importdata('data/L1/pro_90-05_150.txt').data(:,5);
idProIn9 = importdata('data/L1/pro_90-05_150.txt').data(:,1);

% Dissolved Oxygen 88-21
pOxyIn = importdata('data/L1/oxy_88-21_150.txt').data(:,4);
oxyIn = importdata('data/L1/oxy_88-21_150.txt').data(:,5);
idOxyIn = importdata('data/L1/oxy_88-21_150.txt').data(:,1);

% Dissolved Inorganic Carbon 88-21
pDicIn = importdata('data/L1/dic_88-21_150.txt').data(:,4);
dicIn = importdata('data/L1/dic_88-21_150.txt').data(:,5);
idDicIn = importdata('data/L1/dic_88-21_150.txt').data(:,1);

%% Extract Bottle Concentrations within Mixed Layer

% Chlorophyll a
[idChlOut,pChlOut,chlOut] = extractMldVals(idChlIn,pChlIn,chlIn,maxMld);

% Divinyl Chl a
[idDivOut,pDivOut,divOut] = extractMldVals(idDivIn,pDivIn,divIn,maxMld);

% Prochlorococcus: 05-21
[idProOut,pProOut,proOut] = extractMldVals(idProIn,pProIn,proIn,maxMld);

% Prochlorococcus: 90-05
[idProOut9,pProOut9,proOut9] = extractMldVals(idProIn9,pProIn9,proIn9,maxMld);

% Dissolved Oxygen: 88-21
[idOxyOut,pOxyOut,oxyOut] = extractMldVals(idOxyIn,pOxyIn,oxyIn,maxMld);

% Dissolved Inorganic Carbon: 88-21
[idDicOut,pDicOut,dicOut] = extractMldVals(idDicIn,pDicIn,dicIn,maxMld);

%% Visualise ML Extraction

% figure; % Chlorophyll a
% scatter(chlOut,pChlOut);
% grid on;
% set(gca,'YDir','reverse'); ylabel('Pressure [dbar]');
% xlabel('Chl a [ng/l]');
% title('Chl a: 88-21');

% figure; % Divinyl Chl a
% scatter(divOut,pDivOut);
% grid on;
% set(gca,'YDir','reverse'); ylabel('Pressure [dbar]');
% xlabel('Divinyl Chl a [ng/l]');
% title('Divinyl Chl a: 88-21');

% figure; % Prochlorococcus: 05-21
% scatter(proOut,pProOut);
% grid on;
% set(gca,'YDir','reverse'); ylabel('Pressure [dbar]');
% xlabel('$ \textrm{Prochlorococcus [1e5 ml}^{-1}]$','Interpreter','latex');
% title('Prochlorococcus: 05-21');

% figure; % Prochlorococcus: 90-05
% scatter(proOut9,pProOut9);
% grid on;
% set(gca,'YDir','reverse'); ylabel('Pressure [dbar]');
% xlabel('$ \textrm{Prochlorococcus [1e5 ml}^{-1}]$','Interpreter','latex');
% title('Prochlorococcus: 90-05');

% figure; % Oxygen: 88-21
% scatter(oxyOut,pOxyOut);
% grid on;
% set(gca,'YDir','reverse'); ylabel('Pressure [dbar]');
% xlabel('$Oxygen [\mu mol kg}^}{-1}]$','Interpreter','latex');
% title('Oxygen: 88-21');

%% Clean and bin ML extraction

% Chlorophyll a
[~,pChlOutB10,chlOutB,~,~] = cleanAndBin(pChlOut,chlOut,idChlOut');

% Divinyl Chl a
[~,pDivOutB10,divOutB,~,~] = cleanAndBin(pDivOut,divOut,idDivOut');

% Prochlorococcus: 05-21
[~,pProOutB10,proOutB,~,~] = cleanAndBin(pProOut,proOut,idProOut');

% Prochlorococcus: 90-05
[~,pProOutB109,proOutB9,~,~] = cleanAndBin(pProOut9,proOut9,idProOut9');

% Dissolved Oxygen: 88-21
[~,pOxyOutB10,oxyOutB,~,~] = cleanAndBin(pOxyOut,oxyOut,idOxyOut');

% Dissolved Inorganic Carbon: 88-21
[~,pDicOutB10,dicOutB,~,~] = cleanAndBin(pDicOut,dicOut,idDicOut');
%% Find KS p-values, skewness, and kurtosis for ML extraction

% Chlorophyll a
[ksChl,obsChl,pChlKs,chlSk,chlKu,~,~,~] = ksOfBinnedCon(chlOutB,pChlOutB10,10,89);

% Divinyl Chl a
[ksDiv,obsDiv,pDivKs,divSk,divKu,~,~,~] = ksOfBinnedCon(divOutB,pDivOutB10,10,84);

% Prochlorococcus: 05-21
[ksPro,obsPro,pProKs,proSk,proKu,~,~,~] = ksOfBinnedCon(proOutB,pProOutB10,10);

% Prochlorococcus: 90-05
[ksPro9,obsPro9,pProKs9,proSk9,proKu9,~,~,~] = ksOfBinnedCon(proOutB9,pProOutB109,10);

% Dissolved Oxygen: 88-21
[ksOxy,obsOxy,pOxyKs,oxySk,oxyKu,~,~,~] = ksOfBinnedCon(oxyOutB,pOxyOutB10,10,89);

% Dissolved Inorganic Carbon: 88-21
[ksDic,obsDic,pDicKs,dicSk,dicKu,~,~,~] = ksOfBinnedCon(dicOutB,pDicOutB10,10,96);

%% Visualise KS p-values, skewness, and kurtosis

ax1 = figure; % Chlorophyll a
plotKs(pChlKs,ksChl,obsChl,chlSk,chlKu,0.5,20.5,true,89);
sgtitle('[Chl a] 88-21: Mixed Layer');
exportgraphics(ax1,'figures/L1/ks_chla150.png'); clear ax1;

ax2 = figure; % Divinyl Chl a
plotKs(pDivKs,ksDiv,obsDiv,divSk,divKu,0.5,20.5,true,84);
sgtitle('[Divinyl Chl a] 88-21: Mixed Layer');
exportgraphics(ax2,'figures/L1/ks_divi150.png'); clear ax2;

ax3 = figure; % Prochlorococcus: 05-21
plotKs(pProKs,ksPro,obsPro,proSk,proKu,0.5,20.5,true);
sgtitle('Prochlorococcus 05-21: Mixed Layer');
exportgraphics(ax3,'figures/L1/ks_pro150.png'); clear ax3;

ax4 = figure; % Prochlorococcus: 90-05
plotKs(pProKs9,ksPro9,obsPro9,proSk9,proKu9,0.5,20.5,true);
sgtitle('Prochlorococcus 90-05: Mixed Layer');
exportgraphics(ax4,'figures/L1/ks_pro1509.png'); clear ax4;

ax5 = figure; % Dissolved Oxygen: 88-21
plotKs(pOxyKs,ksOxy,obsOxy,oxySk,oxyKu,0.5,20.5,true,89);
sgtitle('Dissolved Oxygen 88-21: Mixed Layer');
exportgraphics(ax5,'figures/L1/ks_oxy150.png'); clear ax5;

ax6 = figure; % Dissolved Inorganic Carbon: 88-21
plotKs(pDicKs,ksDic,obsDic,dicSk,dicKu,0.5,20.5,true,96);
sgtitle('Dissolved Inorganic Carbon 88-21: Mixed Layer');
exportgraphics(ax6,'figures/L1/ks_dic150.png'); clear ax6;

%% Save new data

save mldVals.mat maxMld;
% idChlOut pChlOut idDivOut pDivOut idProOut pProOut;