clear; clc; close all;
addpath("baroneRoutines\"); 
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
set(0,'defaultAxesFontSize',10);

%% 1. Create DCM Array

crn = nan(329,1);
cast = nan(329,31);
pcm = nan(329,31);
tdate = nan(329,31);

ctd = load('datafiles\ctd_iso_ALL.mat').ctd;
for i = 1:329
    tmp = ctd(i).cruise;
    if isempty(tmp)
        tmp = nan;
    end
    crn(i) = tmp;
    cast(i,1:length(ctd(i).cast)) = ctd(i).cast;
    pcm(i,1:length(ctd(i).pcm)) = ctd(i).pcm;
    tdate(i,1:length(ctd(i).date)) = ctd(i).date;
end

% repeat crn
tcrn = repmat(crn',31,1)';

% Combine Above in One Array
tcast = reshape(cast',[],1);
ttcrn = reshape(tcrn',[],1);
tpcm = reshape(pcm',[],1);
tMeanDcm = mean(pcm,2,"omitnan");
tStdDcm = std(pcm,[],2,"omitnan");
tStartCruise = datetime(tdate(:,1),'ConvertFrom','datenum');
ttdate = reshape(tdate',[],1);
tttdate = datetime((ttdate)','ConvertFrom','datenum');
dcmArray = [ttcrn tcast tpcm];
clear i tmp cast crn ctd pcm tcrn tcast ttcrn tpcm ;
% tdate ttdate;

figure
plot(tttdate,dcmArray(:,3));
xlabel('Time'); ylabel('p_{DCM} (dbar)'); set(gca,'YDir','reverse');

figure
errorbar(tStartCruise,tMeanDcm,tStdDcm,'k.'); set(gca,'YDir','reverse');

%%% functionalising from here...
%% 2. Create Bottle Array

% ks = func(botid,pressure,concentration,dcmArray)

bottleId = num2str(importdata('data/hotbot-88_21.txt').data(:,1));
tP = importdata('data/hotbot-88_21.txt').data(:,4);
chla = importdata('data/hotbot-88_21.txt').data(:,5);

id_hplc = num2str(importdata('data\hplcChla_88-21_200.txt').data(:,1));
p_hplc = importdata('data\hplcChla_88-21_200.txt').data(:,4);
hplc = importdata('data\hplcChla_88-21_200.txt').data(:,5);

id_cmo = num2str(importdata('data\chlaMonovinyl_88-21.txt').data(:,1));
p_cmo = importdata('data\chlaMonovinyl_88-21.txt').data(:,4);
cmo = importdata('data\chlaMonovinyl_88-21.txt').data(:,5);

id_cdi = num2str(importdata('data\chlaDivinyl_88-21.txt').data(:,1));
p_cdi = importdata('data\chlaDivinyl_88-21.txt').data(:,4);
cdi = importdata('data\chlaDivinyl_88-21.txt').data(:,5);

id_car = num2str(importdata('data\parC_89-21.txt').data(:,1));
p_car = importdata('data\parC_89-21.txt').data(:,4);
car = importdata('data\parC_89-21.txt').data(:,5);

id_nit = num2str(importdata('data\parN_89-21.txt').data(:,1));
p_nit = importdata('data\parN_89-21.txt').data(:,4);
nit = importdata('data\parN_89-21.txt').data(:,5);

id_pho = num2str(importdata('data\parP_88-21.txt').data(:,1));
p_pho = importdata('data\parP_88-21.txt').data(:,4);
pho = importdata('data\parP_88-21.txt').data(:,5);

id_atp = num2str(importdata('data\atp_88-21.txt').data(:,1));
p_atp = importdata('data\atp_88-21.txt').data(:,4);
atp = importdata('data\atp_88-21.txt').data(:,5);

id_het = num2str(importdata('data\hetBac_90-21.txt').data(:,1));
p_het = importdata('data\hetBac_90-21.txt').data(:,4);
het = importdata('data\hetBac_90-21.txt').data(:,5);

id_pro = num2str(importdata('data\prochl_05-21.txt').data(:,1));
p_pro = importdata('data\prochl_05-21.txt').data(:,4);
pro = importdata('data\prochl_05-21.txt').data(:,5);

id_syn = num2str(importdata('data\synech_05-21.txt').data(:,1));
p_syn = importdata('data\synech_05-21.txt').data(:,4);
syn = importdata('data\synech_05-21.txt').data(:,5);

id_pic = num2str(importdata('data\picoeu_05-21.txt').data(:,1));
p_pic = importdata('data\picoeu_05-21.txt').data(:,4);
pic = importdata('data\picoeu_05-21.txt').data(:,5);

% [ks,obsPerBin,trange] = funShit3(bottleId,tP,chla,dcmArray);        % 88-21
% [tks07,tobsPerBin07,ttrange07] = funShit3(bottleId(1:5890,:),tP(1:5890),chla(1:5890),dcmArray(1:5890,:)); 
% 

%% old -comment
% tCRN = str2num(bottleId(:,1:3)); %str2double doesn't work
% tCast = str2num(bottleId(:,6:8));
% tCast(tCast==100) = nan;

% bottleArray = [tCRN tCast tP];
% clear bottleId tCRN tCast tP;

%% Merge Arrays -comment

% find unique crn/cast combos and remove NaN cast indicators because it is
% unclear to which pcm they apply
% t = rmmissing(unique(bottleArray(:,1:2),"rows"));

% dcmArrayRowNo = []; x = 1;
% for i = 1:length(dcmArray(:,1))
%     if dcmArray(i,1:2) == t(x,1:2)
%         dcmArrayRowNo = [dcmArrayRowNo i];
%         x = x + 1;
%     end
% end

% % save when bottleArray changes
% tid = [];
% for i = 2:9947
%     if bottleArray(i,1) > bottleArray(i-1,1) || bottleArray(i,2) > bottleArray(i-1,2)
%         tid = [tid i];
%     end
% end

% tPcm = nan(9947,1);
% tPcm(1:tid(1)-1) = dcmArray(dcmArrayRowNo(1),3);
% for i = 2:669
%     tPcm(tid(i):tid(i+1)-1) = dcmArray(dcmArrayRowNo(i),3);
% end
% 
% bottleArray = [bottleArray tPcm];

% tPLagrangian = nan(9947,1);
% tPLagrangian = bottleArray(:,3) - bottleArray(:,4);
% bottleArray = [bottleArray tPLagrangian];

% clear i x t dcmArrayRowNo tid tPcm tPLagrangian;

%% Bin pressure -comment
% pB10 = round(bottleArray(:,5),-1);
% bottleArray = [bottleArray pB10];

%% Load chla -comment

% bottleArray = [bottleArray chla];
% tmin = min(bottleArray(:,6));
% tmax = max(bottleArray(:,6));
% trange = tmin:10:tmax;
% clear pB10 chla tmax tmin;
%% KS -comment

% chl = bottleArray(:,7);
% pB = bottleArray(:,6);
% ks = nan(5,41);
% obsPerBin = nan(1,41);
% for i = 1:length(trange)
%     tmp = chl(pB==trange(i));
%     tmp(tmp<=0) = nan;
%     tmp(isnan(tmp)) = [];
%     obsPerBin(i) = length(tmp);
%     if length(tmp) > 3
%         disp(i);
%         [~,ks(:,i),~] = statsplot2(tmp,'noplot');
%     end
% end
% clear tmp;

% chla (regular method)
[bottleArray,trange,chl,pB,ks,obsPerBin] = ksOfLagrangian(bottleId,tP,dcmArray,chla,669);
[bottleArray2,trange07,chl2,pB2,ks07,obsPerBin07] = ksOfLagrangian(bottleId(1:5890,:),tP(1:5890),dcmArray(1:5890,:),chla(1:5890),407);

% chla (HPLC method)
[arrHplc,trHplc,hplcOut,pbHplc,ksHplc,obsPerBinHplc] = ksOfLagrangian(id_hplc,p_hplc,dcmArray,hplc,317);

% HPLC Chlorophyll a (Monovinyl)
[arrCmo,trCmo,cmoOut,pbCmo,ksCmo,obsPerBinCmo] = ksOfLagrangian(id_cmo,p_cmo,dcmArray,cmo,271);

% HPLC Chlorophyll a (Divinyl)
[arrCdi,trCdi,cdiOut,pbCdi,ksCdi,obsPerBinCdi] = ksOfLagrangian(id_cdi,p_cdi,dcmArray,cdi,271);

% Particulate Carbon
[arrCar,trCar,carOut,pbCar,ksCar,obsPerBinCar] = ksOfLagrangian(id_car,p_car,dcmArray,car,332);

% Particulate Nitrogen
[arrNit,trNit,nitOut,pbNit,ksNit,obsPerBinNit] = ksOfLagrangian(id_nit,p_nit,dcmArray,nit,333);

% Particulate Phosphorus
[arrPho,trPho,phoOut,pbPho,ksPho,obsPerBinPho] = ksOfLagrangian(id_pho,p_pho,dcmArray,pho,331);

% ATP - not working, error line 22 in func
% [arrAtp,trAtp,atpOut,pbAtp,ksAtp,obsPerBinAtp] = ksOfLagrangian(id_atp,p_atp,dcmArray,atp,318);

% Heterotrophic Bacteria
[arrHet,trHet,hetOut,pbHet,ksHet,obsPerBinHet] = ksOfLagrangian(id_het,p_het,dcmArray,het,245);

% Prochlorococcus
[arrPro,trPro,proOut,pbPro,ksPro,obsPerBinPro] = ksOfLagrangian(id_pro,p_pro,dcmArray,pro,147);

% Synechococcus
[arrSyn,trSyn,synOut,pbSyn,ksSyn,obsPerBinSyn] = ksOfLagrangian(id_syn,p_syn,dcmArray,syn,147);

% Picoeukaryotes
[arrPic,trPic,picOut,pbPic,ksPic,obsPerBinPic] = ksOfLagrangian(id_pic,p_pic,dcmArray,pic,147);

%% Chlorophyll a (Regular Method): Lagrangian 10 dbar

ax1 = figure;
subplot(1,2,1)
barh(obsPerBin);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([10 35]);
set(gca,"YTick",1:1:38,"YTickLabel",-240:10:130);
title('No. of Observations');

subplot(1,2,2)
plot(ks(1,:),trange,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ks(2,:),trange,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks(3,:),trange,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks(4,:),trange,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
xlim([0 0.8]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle chl-a [10 dbar bins] (88-21)');
exportgraphics(ax1,'figures/ks_bottleLagrangian10db.png');
clear ax1;

%% Chlorophyll a (Regular Method): Lagrangian 10 dbar: 88-07 -comment

ax2 = figure;
subplot(1,2,1)
barh(obsPerBin07);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([10 35]);
set(gca,"YTick",1:1:length(trange07),"YTickLabel",trange07);
title('No. of Observations');

subplot(1,2,2)
plot(ks07(1,:),trange07,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ks07(2,:),trange07,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks07(3,:),trange07,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks07(4,:),trange07,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle chl-a [10 dbar bins] (88-07)');
exportgraphics(ax2,'figures/ks_bottleLagrangian10db07.png');
clear ax2;

%% HPLC Chlorophyll a: Lagrangian 10 dbar

ax3 = figure;
subplot(1,2,1)
barh(obsPerBinHplc);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([1 26]);
set(gca,"YTick",1:1:length(trHplc),"YTickLabel",trHplc);
title('No. of Observations');

subplot(1,2,2)
plot(ksHplc(1,:),trHplc,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksHplc(2,:),trHplc,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksHplc(3,:),trHplc,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksHplc(4,:),trHplc,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('HPLC Method: [chl a] (Lagrangian, 10 dbar bin)');
exportgraphics(ax3,'figures/ks_HplcLag.png');
clear ax3;

%% HPLC Monovinyl Chlorophyll a: Lagrangian 10 dbar

ax4 = figure;
subplot(1,2,1)
barh(obsPerBinCmo);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([1 26]);
set(gca,"YTick",1:1:length(trCmo),"YTickLabel",trCmo);
title('No. of Observations');

subplot(1,2,2)
plot(ksCmo(1,:),trCmo,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksCmo(2,:),trCmo,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksCmo(3,:),trCmo,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksCmo(4,:),trCmo,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('HPLC Monovinyl Chlorophyll a (Lagrangian, 10 dbar bin)');
exportgraphics(ax4,'figures/ks_CmoLag.png');
clear ax4;

%% HPLC Divinyl Chlorophyll a: Lagrangian 10 dbar

ax5 = figure;
subplot(1,2,1)
barh(obsPerBinCdi);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([1 26]);
set(gca,"YTick",1:1:length(trCdi),"YTickLabel",trCdi);
title('No. of Observations');

subplot(1,2,2)
plot(ksCdi(1,:),trCdi,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksCdi(2,:),trCdi,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksCdi(3,:),trCdi,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksCdi(4,:),trCdi,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('HPLC Divinyl Chlorophyll a (Lagrangian, 10 dbar bin)');
exportgraphics(ax5,'figures/ks_CdiLag.png');
clear ax5;

%% Particulate Carbon: Lagrangian 10 dbar

ax6 = figure;
subplot(1,2,1)
barh(obsPerBinCar);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([3 28]);
set(gca,"YTick",1:1:length(trCar),"YTickLabel",trCar);
title('No. of Observations');

subplot(1,2,2)
plot(ksCar(1,:),trCar,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksCar(2,:),trCar,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksCar(3,:),trCar,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksCar(4,:),trCar,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Particulate Carbon (Lagrangian, 10 dbar bin)');
exportgraphics(ax6,'figures/ks_CarLag.png');
clear ax6;

%% Particulate Nitrogen: Lagrangian 10 dbar

ax7 = figure;
subplot(1,2,1)
barh(obsPerBinNit);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([3 28]);
set(gca,"YTick",1:1:length(trNit),"YTickLabel",trNit);
title('No. of Observations');

subplot(1,2,2)
plot(ksNit(1,:),trNit,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksNit(2,:),trNit,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksNit(3,:),trNit,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksNit(4,:),trNit,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Particulate Nitrogen (Lagrangian, 10 dbar bin)');
exportgraphics(ax7,'figures/ks_NitLag.png');
clear ax7;

%% Particulate Phosphorus: Lagrangian 10 dbar

ax8 = figure;
subplot(1,2,1)
barh(obsPerBinPho);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([3 28]);
set(gca,"YTick",1:1:length(trPho),"YTickLabel",trPho);
title('No. of Observations');

subplot(1,2,2)
plot(ksPho(1,:),trPho,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPho(2,:),trPho,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPho(3,:),trPho,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPho(4,:),trPho,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Particulate Phosphorus (Lagrangian, 10 dbar bin)');
exportgraphics(ax8,'figures/ks_PhoLag.png');
clear ax8;

%% Heterotrophic Bacteria: Lagrangian 10 dbar

ax9 = figure;
subplot(1,2,1)
barh(obsPerBinHet);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([3 28]);
set(gca,"YTick",1:1:length(trHet),"YTickLabel",trPho);
title('No. of Observations');

subplot(1,2,2)
plot(ksHet(1,:),trHet,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksHet(2,:),trHet,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksHet(3,:),trHet,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksHet(4,:),trHet,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Heterotrophic Bacteria (Lagrangian, 10 dbar bin)');
exportgraphics(ax9,'figures/ks_HetLag.png');
clear ax9;

%% Prochlorococcus: Lagrangian 10 dbar

ax10 = figure;
subplot(1,2,1)
barh(obsPerBinPro);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
% ylim([3 28]);
set(gca,"YTick",1:1:length(trPro),"YTickLabel",trPro);
title('No. of Observations');

subplot(1,2,2)
plot(ksPro(1,:),trPro,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPro(2,:),trPro,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPro(3,:),trPro,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPro(4,:),trPro,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Prochlorococcus (Lagrangian, 10 dbar bin)');
exportgraphics(ax10,'figures/ks_ProLag.png');
clear ax10;

%% Synechococcus: Lagrangian 10 dbar

ax11 = figure;
subplot(1,2,1)
barh(obsPerBinSyn);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
% ylim([3 28]);
set(gca,"YTick",1:1:length(trSyn),"YTickLabel",trSyn);
title('No. of Observations');

subplot(1,2,2)
plot(ksSyn(1,:),trSyn,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksSyn(2,:),trSyn,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksSyn(3,:),trSyn,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksSyn(4,:),trSyn,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Synechococcus (Lagrangian, 10 dbar bin)');
exportgraphics(ax11,'figures/ks_SynLag.png');
clear ax11;

%% Picoeukaryotes: Lagrangian 10 dbar

ax12 = figure;
subplot(1,2,1)
barh(obsPerBinPic);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
% ylim([3 28]);
set(gca,"YTick",1:1:length(trPic),"YTickLabel",trPic);
title('No. of Observations');

subplot(1,2,2)
plot(ksPic(1,:),trPic,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksPic(2,:),trPic,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksPic(3,:),trPic,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksPic(4,:),trPic,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 100]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Picoeukaryotes (Lagrangian, 10 dbar bin)');
exportgraphics(ax12,'figures/ks_PicLag.png');
clear ax12;
%% DELETE ME later
% The following code was to calculate the Lagrangian KS for chla (regular
% method) before I functionalised it properly.

% chl07 = chl(1:5890);
% pB07 = pB(1:5890);
% ks07 = nan(5,41);
% obsPerBin07 = nan(1,41);
% 
% for i = 1:length(trange)
%     tmp1 = chl07(pB07==trange(i));
%     tmp1(tmp1<=0) = nan;
%     tmp1(isnan(tmp1)) = [];
%     obsPerBin07(i) = length(tmp1);
%     if length(tmp1) > 3
%         disp(i);
%         [~,ks07(:,i),~] = statsplot2(tmp1,'noplot');
%     end
% end
% clear tmp1;
% 
% for i = 1:41
%     if obsPerBin07(i) < 100
%         ks07(:,i) = nan;
%     end
% end
