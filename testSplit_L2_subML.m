clear; clc; close all; addpath("baroneRoutines\");

p_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,4);
chl_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,5);
id_hplc = num2str(importdata('data/HPLC_chla_88-21.txt').data(:,1));

% Load CTD Salinity and Temperature: save where readings coincide with chl
tmpT = importdata('data\HPLC_chlaTS_88-21.txt').data(:,5);
tmpSp = importdata('data\HPLC_chlaTS_88-21.txt').data(:,6);
tmpChl = importdata('data\HPLC_chlaTS_88-21.txt').data(:,7);
tmpChl(tmpChl==-9) = nan;
t_hplc = tmpT(~isnan(tmpChl));
sp_hplc = tmpSp(~isnan(tmpChl));

clear tmpT tmpSp tmpChl;

%%

% Only look at values below the maximum MLD of a cruise
% We don't have a convenient measurement of MLD in 'isopycnal' coordinates.
pMaxMld = load('testPMld.mat').maxMldPerCruise;

bottleCRN = str2num(id_hplc(:,1:3));

%%

% these new values are pressures etc that are below MLD

L = 3550;

tmpP_subML = nan(L,1);
tmpCRN_subML = nan(L,1);
tmpChl_subML = nan(L,1);
% tmpId_subML = nan(L,1);
tmpCrnStr_subML = nan(L,1);
tmpT = nan(L,1); tmpS = nan(L,1);
botID = [];

for i = 1:L
    %bottleCRN(i)
    % with MAX MLD per cruise
    tmpMld = pMaxMld(bottleCRN(i));
    if p_hplc(i) > tmpMld
        tmpP_subML(i) = p_hplc(i);
        tmpCRN_subML(i) = bottleCRN(i);
        tmpChl_subML(i) = chl_hplc(i);
        tmpStr = id_hplc(i,:);
        tmpT(i) = t_hplc(i);
        tmpS(i) = sp_hplc(i);
        botID = [botID;tmpStr];
    end
end

% now what do I do...
pSubML = tmpP_subML(~isnan(tmpP_subML));
crnSubML = tmpCRN_subML(~isnan(tmpCRN_subML));
chlSubML = tmpChl_subML(~isnan(tmpChl_subML));
% idSubML = tmpId_subML(~isnan(tmpId_subML));
tSubML = tmpT(~isnan(tmpT));
sSubML = tmpS(~isnan(tmpS));

%% Calculate Density for above

SA_subML = gsw_SA_from_SP(sSubML,pSubML,158,22.75);
CT_subML = gsw_CT_from_t(SA_subML,tSubML,pSubML);

sigma0_subML = gsw_sigma0(SA_subML,CT_subML);

sigma0_subML_ROUND = round(sigma0_subML,5);

%%
figure;
histogram(sigma0_subML_ROUND);

%%
figure
scatter(chlSubML,pSubML,'Marker','.');
set(gca,'YDir','reverse'); ylim([0 200]);
xlabel('Chl a (HPLC) [ng/l]'); ylabel('Pressure [dbar]');
title('Sub-ML chl a','Interpreter','latex');

%%


dcm = load("dcm.mat").dcm;

%% Calculate KS p-value for Chl-a (HPLC)

% Threshold = 100
[p100,ks100,obs100,sk100,ku100,~,~,~] = ksOfLagrangian(botID,pSubML,dcm,chlSubML,159,100);
% Threshold = 80
[p80,ks80,obs80,sk80,ku80,~,~,~] = ksOfLagrangian(botID,pSubML,dcm,chlSubML,159,77);
% Threshold = 48
[p45,ks45,obs45,sk45,ku45,~,~,~] = ksOfLagrangian(botID,pSubML,dcm,chlSubML,159,48);


%% Plot KS p-value for Chl-a (HPLC)

% Threshold = 100
ax1 = figure;
plotKs(p100,ks100,obs100,sk100,ku100,1,23,false,100,[-120 100]);
sgtitle('KS Test: Chl-a (HPLC), T = 100, 88-21');
exportgraphics(ax1,'figures/ks_botLagSubML_t100.png'); clear ax1;

% Threshold = 77
ax2 = figure;
plotKs(p80,ks80,obs80,sk80,ku80,1,23,false,77,[-120 100]);
sgtitle('KS Test: Chl-a (HPLC), T = 80, 88-21');
exportgraphics(ax2,'figures/ks_botLagSubML_t80.png'); clear ax2;

% Threshold = 48
ax3 = figure;
plotKs(p45,ks45,obs45,sk45,ku45,1,23,false,48,[-120 100]);
sgtitle('KS Test: Chl-a (HPLC), T = 45, 88-21');
exportgraphics(ax3,'figures/ks_botLagSubML_t45.png'); clear ax3;



%% density check for cruise no. 2

% test = importdata('data/crn2_chlTS.txt').data(:,4);
testP = importdata('data/crn2_chlTS.txt').data(:,4);
testT = importdata('data\crn2_chlTS.txt').data(:,5);
testS_b = importdata('data\crn2_chlTS.txt').data(:,7);
testS_c = importdata('data/crn2_chlTS.txt').data(:,6);
testChl = importdata('data/crn2_chlTS.txt').data(:,8);
testId = num2str(importdata('data/crn2_chlTS.txt').data(:,1));

cast = str2num(testId(:,5:6));

SA = gsw_SA_from_SP([testS_c],testP,158,22.75);
SA(SA==0) = nan;
CT = gsw_CT_from_t(SA,testT,testP);

siggyBoi = gsw_sigma0(SA,CT);

newVec = [cast(35:43) testP(35:43) siggyBoi(35:43)];



%%
sigCtd = [23.6888 24.1807 24.6807 25.0348];

ax4 = figure;
yyaxis left; scatter(testChl(40:43),sigCtd); ylabel('$\sigma_0$ (in-situ)','Interpreter','latex'); hold on
yyaxis right; scatter(testChl(40:43),siggyBoi(40:43)); ylabel('$\sigma_0$: TEOS-10','Interpreter','latex'); hold off
xlabel('HPLC Chl a [ng/l]'); title('Density of chl a measurements: Cruise No. 2 (HOT, 1988)');
exportgraphics(ax4,'figures/crn2_sigma0_comparison.png'); clear ax4;

%% density: all cruises

% test = importdata('data/crn2_chlTS.txt').data(:,4);
testP = importdata('data/chlTS.txt').data(:,4);
testT = importdata('data/chlTS.txt').data(:,5);
% testS_b = importdata('data/chlTS.txt').data(:,7);
testS_c = importdata('data/chlTS.txt').data(:,6);
testChl = importdata('data/chlTS.txt').data(:,7);
testId = num2str(importdata('data/chlTS.txt').data(:,1));
testChl(testChl == -9) = nan;

cast = str2num(testId(:,5:6));

SA = gsw_SA_from_SP(testS_c,testP,158,22.75);
SA(SA==0) = nan;
CT = gsw_CT_from_t(SA,testT,testP);

siggyBoi = gsw_sigma0(SA,CT);

% newVec = [cast(35:43) testP(35:43) siggyBoi(35:43)];
figure;scatter(testChl,siggyBoi); set(gca,'YDir','reverse');

sigBin = discretize(siggyBoi,22.4:0.1:26);

% sigBin2 = discretize(siggyBoi,min(siggyBoi):0.1:max(siggyBoi));
%%
% [ks, obs, depth2, Sk, ku, sd, c95, mu] = ksOfBinnedCon(testChl, sigBin, 10);

ksSig = nan(5,36);
Sk = nan(1,36); Ku = nan(1,36); obs = nan(1,36);

for i = 1:36
    X_i = testChl(sigBin==i);
    X_i(isnan(X_i)) = [];
    if length(X_i) >3
        [~,ksSig(:,i),~,~,~,~] = statsplot2(X_i,'noplot');
        Sk(i) = skewness(X_i);
        Ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
end

% figure;
% plot(ksSig,1:1:36);
% set(gca,'YDir','reverse');

tr = linspace(min(siggyBoi),max(siggyBoi),36);

plotKsSig(tr,ksSig,obs,Sk,Ku);