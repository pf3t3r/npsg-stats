clear; clc; close all; addpath("baroneRoutines\");

%% Load MLD and DCM
pMaxMld = load('testPMld.mat').maxMldPerCruise;
dcm = load("dcm.mat").dcm;  % pDcm and sigmaDcm (all casts, crn 1 - 329)

%% Load parameters to test

% CHL
% Load bottle ID, pressure, chl, and corresponding CTD Sp and T
chlId = num2str(importdata('data/HPLC_chlaTS_88-21.txt').data(:,1));
chlP = importdata('data/HPLC_chlaTS_88-21.txt').data(:,4);
chlChl = importdata('data\HPLC_chlaTS_88-21.txt').data(:,7);
chlT = importdata('data\HPLC_chlaTS_88-21.txt').data(:,5);
chlSp = importdata('data\HPLC_chlaTS_88-21.txt').data(:,6);

% Let's try same for DIVINYL CHL A
divId = num2str(importdata('data/chlaDivinyl_TS_88-21.txt').data(:,1));
divP = importdata('data/chlaDivinyl_TS_88-21.txt').data(:,4);
divDiv = importdata('data/chlaDivinyl_TS_88-21.txt').data(:,7);
divT = importdata('data/chlaDivinyl_TS_88-21.txt').data(:,5);
divSp = importdata('data/chlaDivinyl_TS_88-21.txt').data(:,6);

% Let's try for PROCHLOROCOCCUS 05-21
proId = num2str(importdata('data/proTS_05-21.txt').data(:,1));
proP = importdata('data/proTS_05-21.txt').data(:,4);
proPro = importdata('data/proTS_05-21.txt').data(:,7);
proT = importdata('data/proTS_05-21.txt').data(:,5);
proSp = importdata('data/proTS_05-21.txt').data(:,6);


%% Find where density and concentration measurements coincide

% Chl a
[chlIdOut,sigChl,chlOut,testcrn] = densityConcentrations(chlId, chlP, chlChl, chlT, chlSp, pMaxMld);

% Divinyl Chl a
[divIdOut,sigDiv,divOut] = densityConcentrations(divId, divP, divDiv, divT, divSp, pMaxMld);

% Prochlorococcus
[proIdOut,sigPro,proOut] = densityConcentrations(proId, proP, proPro, proT, proSp, pMaxMld);

%% Apply ksOfIsoLagrangian function

t93 = 93; t67 = 67; t39 = 39; t31 = 31;

% CHLA (HPLC)
[trange,ks,obsPerBin,Sk,Ku,botArr,sigB,Xout] = ksOfIsoLagrangian(chlIdOut,sigChl,dcm,chlOut,t93);
[tr67,ks67,obs67,sk67,ku67,bot67,sigB67,Xout67] = ksOfIsoLagrangian(chlIdOut,sigChl,dcm,chlOut,t67);
[tr39,ks39,obs39,sk39,ku39,bot39,sigB39,Xout39] = ksOfIsoLagrangian(chlIdOut,sigChl,dcm,chlOut,t39);

% DIVINYL CHLA (HPLC)
[trDiv,ksDiv,obsDiv,SkDiv,KuDiv,botDiv,sigBDiv,XoutDiv] = ksOfIsoLagrangian(divIdOut,sigDiv,dcm,divOut,t93);

% Prochlorococcus
[trPro,ksPro,obsPro,SkPro,KuPro,botPro,sigBPro,XoutPro] = ksOfIsoLagrangian(proIdOut,sigPro,dcm,proOut,t31);



%
% axA = figure;
% plotKsSig(obsPerBin,t93,[1 38],[-1.9 1.8],ks,trange,Sk,Ku);
% sgtitle('Chl-a (HPLC) 88-21: sub-ML, isopycnal, DCM-centred','FontSize',10);
% exportgraphics(axA,'figures/L2/ks_chla_func.png'); clear axA;
% %
% % Extract id, p, chl, T, Sp only where there are measurements of chl
% tmpChl(tmpChl==-9) = nan;
% 
% id_hplc = chlId(~isnan(tmpChl),:);
% p_hplc = chlP(~isnan(tmpChl));
% chl_hplc = tmpChl(~isnan(tmpChl));
% t_hplc = tmpT(~isnan(tmpChl));
% sp_hplc = tmpSp(~isnan(tmpChl));
% 
% clear tmpT tmpSp;
%
% % Only look at values below the maximum MLD of a cruise
% % We don't have a convenient measurement of MLD in 'isopycnal' coordinates.
% 
% % Load mixed layer depth 'pMaxMld'. It is defined as the maximum mixed
% % layer depth attained during a cruise.
% pMaxMld = load('testPMld.mat').maxMldPerCruise;
% 
% % Extract cruise number 'crn'
% crn = str2num(id_hplc(:,1:3));
% Extract measurements taken below pMaxMld
% L = 3550; % No. of casts with chl measurements in cruise 1 - 329
% 
% tmpP_subML = nan(L,1);
% tmpCRN_subML = nan(L,1);
% tmpChl_subML = nan(L,1);
% % tmpCrnStr_subML = nan(L,1);
% tmpT = nan(L,1); tmpS = nan(L,1);
% botID = [];
% 
% % To find the measurements taken beneath pMaxMld, we save them to a new
% % file...
% for i = 1:L
%     %crn(i)
%     % with MAX MLD per cruise
%     tmpMld = pMaxMld(crn(i));
%     if p_hplc(i) > tmpMld
%         tmpP_subML(i) = p_hplc(i);
%         tmpCRN_subML(i) = crn(i);
%         tmpChl_subML(i) = chl_hplc(i);
%         tmpStr = id_hplc(i,:);
%         tmpT(i) = t_hplc(i);
%         tmpS(i) = sp_hplc(i);
%         botID = [botID;tmpStr];
%     end
% end
% 
% % ...and remove nan values (which here represent measurements above
% % pMaxMld)
% pSubML = tmpP_subML(~isnan(tmpP_subML));
% % crnSubML = tmpCRN_subML(~isnan(tmpCRN_subML));
% chlSubML = tmpChl_subML(~isnan(tmpChl_subML));
% tSubML = tmpT(~isnan(tmpT));
% sSubML = tmpS(~isnan(tmpS));
% % Calculate the sub-ML Density 'sig0_sM'
% SA_subML = gsw_SA_from_SP(sSubML,pSubML,158,22.75);
% CT_subML = gsw_CT_from_t(SA_subML,tSubML,pSubML);
% 
% sig0_sM = gsw_sigma0(SA_subML,CT_subML);
% % Bin the sub-ML density -> 'sig0_sMr'
% sig0_sMr = round(sig0_sM,3,"significant");
% % Display Histogram of the binned sub-ML density 'sig0_sMr' 
% % Edges = (no. of unique values/bins + missing bins) - 1
% uniH = length(unique(sig0_sMr));
% edgesH = length(min(sig0_sMr):0.1:max(sig0_sMr)) - 1;
% 
% figure;
% histogram(sig0_sMr,edgesH);
% % Remove densest waters -> continuous, binned sub-ML density 'sig0_sMrf'
% This is justified because (1) it will give a continuous 'depth' profile
% of density, and (2) 99.4% of the measurements are contained within the
% continuosu portion anyway.
% sig0_sMrf = sig0_sMr(sig0_sMr<=25.8);
% 
% % Confirm that filtering has been done correctly
% figure; histogram(sig0_sMrf,'BinWidth',0.1,'BinLimits',[min(sig0_sMrf) max(sig0_sMrf)]);
% % Extract Chl-a at the filtered, binned densities found above
% chlaSubMLf = chlSubML(sig0_sMr<=25.8);
% botIDf = botID(sig0_sMr<=25.8,:);             % this is needed for ksOfLagrangian()
% % crnSubMLf = crnSubML(sig0_sMr<=25.8);
% % Load DCM vector
% This shows pDcm and sigmaDcm of each cast for cruises 1 - 329
% dcm = load("dcm.mat").dcm;
% % Apply ksOfIsoLagrangian function
% t93 = 93; t67 = 67; t39 = 39;
% 
% [trange,ks,obsPerBin,Sk,Ku,botArr,sigB,Xout] = ksOfIsoLagrangian(botIDf,sig0_sMrf,dcm,chlaSubMLf,t93);
% Alternative thresholds: 67, 39
% [tr67,ks67,obs67,sk67,ku67,bot67,sigB67,Xout67] = ksOfIsoLagrangian(botIDf,sig0_sMrf,dcm,chlaSubMLf,t67);
% [tr39,ks39,obs39,sk39,ku39,bot39,sigB39,Xout39] = ksOfIsoLagrangian(botIDf,sig0_sMrf,dcm,chlaSubMLf,t39);

%% Chl a: Compare untransformed coordinates (Eulerian) to DCM-centred, Isopycnal coordinates

figure;
yyaxis left
scatter(chlOut,sigChl,'Marker','.'); set(gca,'YDir','reverse');
ylabel('$\sigma_0: \textrm{below ML } [\textrm{kg m}^{-3}]$','Interpreter','latex');
yyaxis right
scatter(botArr(:,7),botArr(:,6),'Marker','.');
set(gca,'YDir','reverse'); xlabel('$\textrm{Chl } \textit{a } [\textrm{ng l}^{-1}$]','Interpreter','latex');
ylabel('$\sigma_0: \textrm{DCM-centred, below ML} [\textrm{kg m}^{-3}]$','Interpreter','latex');
title('Chl a vs. $\sigma_0$',[],'Interpreter','latex');

%% KS p-values: Chl a (Isopycnal, DCM-centred)

% trange(20) = round(trange(20));

ax1 = figure;
plotKsSig(obsPerBin,t93,[1 38],[-1.9 1.8],ks,trange,Sk,Ku);
sgtitle('Chl-a (HPLC) 88-21: sub-ML, isopycnal, DCM-centred','FontSize',10);
exportgraphics(ax1,'figures/L2/ks_chla.png'); clear ax1;

% ax1a = figure;
% plotKsSig(obsPerBin,t67,[1 38],[-1.9 1.8],ks67,tr67,sk67,ku67);
% sgtitle('Chl-a (HPLC) 88-21: sub-ML, isopycnal, DCM-centred (t = 67)','FontSize',10);
% exportgraphics(ax1a,'figures/L2/ks_chla_t67.png'); clear ax1a;
% 
% ax1b = figure;
% plotKsSig(obsPerBin,t39,[1 38],[-1.9 1.8],ks39,tr39,sk39,ku39);
% sgtitle('Chl-a (HPLC) 88-21: sub-ML, isopycnal, DCM-centred (t = 39)','FontSize',10);
% exportgraphics(ax1b,'figures/L2/ks_chla_t39.png'); clear ax1b;

%% KS p-values: Divinyl Chl a (Isopycnal, DCM-centred)

ax2 = figure;
plotKsSig(obsDiv,t93,[1 38],[-1.9 1.8],ksDiv,trDiv,SkDiv,KuDiv);
sgtitle('Divinyl Chl-a (HPLC) 88-21: sub-ML, isopycnal, DCM-centred','FontSize',10);
exportgraphics(ax2,'figures/L2/ks_diviChla.png'); clear ax2;

%% KS p-values: Prochlorococcus (Isopycnal, DCM-centred)

ax3 = figure;
plotKsSig(obsPro,t31,[1 38],[-1.9 1.8],ksPro,trPro,SkPro,KuPro);
sgtitle('Prochlorococcus (05-21): sub-ML, isopycnal, DCM-centred','FontSize',10);
exportgraphics(ax3,'figures/L2/ks_pro.png'); clear ax3;
