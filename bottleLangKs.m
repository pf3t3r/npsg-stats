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
dcm = [ttcrn tcast tpcm];

figure
plot(tttdate,dcm(:,3));
xlabel('Time'); ylabel('p_{DCM} (dbar)'); set(gca,'YDir','reverse');

figure
errorbar(tStartCruise,tMeanDcm,tStdDcm,'k.'); set(gca,'YDir','reverse');

clear i tmp cast crn ctd pcm tcrn tcast ttcrn tpcm tdate ttdate tttdate;

%% 2. Load bottle concentrations (id, pressure, concentration)

idChla = num2str(importdata('data/hotbot-88_21.txt').data(:,1));
pChla = importdata('data/hotbot-88_21.txt').data(:,4);
chla = importdata('data/hotbot-88_21.txt').data(:,5);

idHplc = num2str(importdata('data\hplcChla_88-21_200.txt').data(:,1));
pHplc = importdata('data\hplcChla_88-21_200.txt').data(:,4);
hplc = importdata('data\hplcChla_88-21_200.txt').data(:,5);

idCmo = num2str(importdata('data\chlaMonovinyl_88-21.txt').data(:,1));
pCmo = importdata('data\chlaMonovinyl_88-21.txt').data(:,4);
cmo = importdata('data\chlaMonovinyl_88-21.txt').data(:,5);

idCdi = num2str(importdata('data\chlaDivinyl_88-21.txt').data(:,1));
pCdi = importdata('data\chlaDivinyl_88-21.txt').data(:,4);
cdi = importdata('data\chlaDivinyl_88-21.txt').data(:,5);

idCar = num2str(importdata('data\parC_89-21.txt').data(:,1));
pCar = importdata('data\parC_89-21.txt').data(:,4);
car = importdata('data\parC_89-21.txt').data(:,5);

idNit = num2str(importdata('data\parN_89-21.txt').data(:,1));
pNit = importdata('data\parN_89-21.txt').data(:,4);
nit = importdata('data\parN_89-21.txt').data(:,5);

idPho = num2str(importdata('data\parP_88-21.txt').data(:,1));
pPho = importdata('data\parP_88-21.txt').data(:,4);
pho = importdata('data\parP_88-21.txt').data(:,5);

idAtp = num2str(importdata('data\atp_88-21.txt').data(:,1));
pAtp = importdata('data\atp_88-21.txt').data(:,4);
atp = importdata('data\atp_88-21.txt').data(:,5);

idHet = num2str(importdata('data\hetBac_90-21.txt').data(:,1));
pHet = importdata('data\hetBac_90-21.txt').data(:,4);
het = importdata('data\hetBac_90-21.txt').data(:,5);

idPro = num2str(importdata('data\prochl_05-21.txt').data(:,1));
pPro = importdata('data\prochl_05-21.txt').data(:,4);
pro = importdata('data\prochl_05-21.txt').data(:,5);

idSyn = num2str(importdata('data\synech_05-21.txt').data(:,1));
pSyn = importdata('data\synech_05-21.txt').data(:,4);
syn = importdata('data\synech_05-21.txt').data(:,5);

idPic = num2str(importdata('data\picoeu_05-21.txt').data(:,1));
pPic = importdata('data\picoeu_05-21.txt').data(:,4);
pic = importdata('data\picoeu_05-21.txt').data(:,5);

%% 3. Apply Kolmogorov Smirnov Test to Bottle Concentrations

% Default Threshold (=100)

[trChla,ksChla,obsChla] = ksOfLagrangian(idChla,pChla,dcm,chla,669);                      % chla (regular method)
b = 1:5890;
[trChla07,ksChla07,obsChla07] = ksOfLagrangian(idChla(b,:),pChla(b),dcm(b,:),chla(b),407);    % chla (regular method), 1988-2007
[trHplc,ksHplc,obsHplc] = ksOfLagrangian(idHplc,pHplc,dcm,hplc,317);                        % chla (HPLC method)
[trCmo,ksCmo,obsCmo] = ksOfLagrangian(idCmo,pCmo,dcm,cmo,271);                              % HPLC Chlorophyll a (Monovinyl)
[trCdi,ksCdi,obsCdi] = ksOfLagrangian(idCdi,pCdi,dcm,cdi,271);                              % HPLC Chlorophyll a (Divinyl)
[trCar,ksCar,obsCar] = ksOfLagrangian(idCar,pCar,dcm,car,332);                              % Particulate Carbon
[trNit,ksNit,obsNit] = ksOfLagrangian(idNit,pNit,dcm,nit,333);                              % Particulate Nitrogen
[trPho,ksPho,obsPho] = ksOfLagrangian(idPho,pPho,dcm,pho,331);                              % Particulate Phosphorus

% ATP - not working, error line 22 in func
% [arrAtp,trAtp,atpOut,pbAtp,ksAtp,obsPerBinAtp] = ksOfLagrangian(idAtp,pAtp,dcmArray,atp,318);

[trHet,ksHet,obsHet] = ksOfLagrangian(idHet,pHet,dcm,het,245);                              % Heterotrophic Bacteria
[trPro,ksPro,obsPro] = ksOfLagrangian(idPro,pPro,dcm,pro,147);                              % Prochlorococcus
[trSyn,ksSyn,obsSyn] = ksOfLagrangian(idSyn,pSyn,dcm,syn,147);                              % Synechococcus
[trPic,ksPic,obsPic] = ksOfLagrangian(idPic,pPic,dcm,pic,147);                              % Picoeukaryotes

% Lower Threshold (=50)
[trPro50,ksPro50,obsPro50] = ksOfLagrangian(idPro,pPro,dcm,pro,147,50);                              % Prochlorococcus
[trSyn50,ksSyn50,obsSyn50] = ksOfLagrangian(idSyn,pSyn,dcm,syn,147,50);                              % Synechococcus
[trPic50,ksPic50,obsPic50] = ksOfLagrangian(idPic,pPic,dcm,pic,147,50);                              % Picoeukaryotes

%% 4. Plot the KS Statistics for each Bottle Concentration Time Series
% All of the following results are DCM-centred (Lagrangian), binned to 10
% dbar depths, and comprise the full HOT time series from 1988 up to 2021
% (2022- was not ready at time of download) unless otherwise specified.

%% 4.a.i. Chlorophyll a (Regular Method)

ax1 = figure;
plotKs(trChla,ksChla,obsChla,10,35,false);
sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle chl-a [10 dbar bins] (88-21)');
exportgraphics(ax1,'figures/ks_bottleLagrangian10db.png'); clear ax1;

%% 4.a.ii. Chlorophyll a (Regular Method): 88-07

ax2 = figure;
plotKs(trChla07,ksChla07,obsChla07,10,35,false);
sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle chl-a [10 dbar bins] (88-07)');
exportgraphics(ax2,'figures/ks_bottleLagrangian10db07.png'); clear ax2;

%% 4.b. HPLC Chlorophyll a

ax3 = figure;
plotKs(trHplc,ksHplc,obsHplc,1,26,false);
sgtitle('HPLC Method: [chl a] (Lagrangian, 10 dbar bin)');
exportgraphics(ax3,'figures/ks_HplcLag.png'); clear ax3;

%% 4.c. HPLC Monovinyl Chlorophyll a

ax4 = figure;
plotKs(trCmo,ksCmo,obsCmo,1,26,false);
sgtitle('HPLC Monovinyl Chlorophyll a (Lagrangian, 10 dbar bin)');
exportgraphics(ax4,'figures/ks_CmoLag.png'); clear ax4;

%% 4.d. HPLC Divinyl Chlorophyll a

ax5 = figure;
plotKs(trCdi,ksCdi,obsCdi,1,26,false);
sgtitle('HPLC Divinyl Chlorophyll a (Lagrangian, 10 dbar bin)');
exportgraphics(ax5,'figures/ks_CdiLag.png'); clear ax5;

%% 4.e. Particulate Carbon

ax6 = figure;
plotKs(trCar,ksCar,obsCar,3,28,false);
sgtitle('Particulate Carbon (Lagrangian, 10 dbar bin)');
exportgraphics(ax6,'figures/ks_CarLag.png'); clear ax6;

%% 4.f. Particulate Nitrogen

ax7 = figure;
plotKs(trNit,ksNit,obsNit,3,28,false);
sgtitle('Particulate Nitrogen (Lagrangian, 10 dbar bin)');
exportgraphics(ax7,'figures/ks_NitLag.png'); clear ax7;

%% 4.g. Particulate Phosphorus

ax8 = figure;
plotKs(trPho,ksPho,obsPho,3,28,false);
sgtitle('Particulate Phosphorus (Lagrangian, 10 dbar bin)');
exportgraphics(ax8,'figures/ks_PhoLag.png'); clear ax8;

%% 4.h. Heterotrophic Bacteria

ax9 = figure;
plotKs(trHet,ksHet,obsHet,1,26,false);
sgtitle('Heterotrophic Bacteria (Lagrangian, 10 dbar bin)');
exportgraphics(ax9,'figures/ks_HetLag.png'); clear ax9;

%% 4.i.i Prochlorococcus

ax10 = figure;
plotKs(trPro,ksPro,obsPro,1,26,false);
sgtitle('Prochlorococcus (Lagrangian, 10 dbar bin)');
exportgraphics(ax10,'figures/ks_ProLag.png'); clear ax10;

%% 4.i.ii Prochlorococcus (threshold = 50)

ax10a = figure;
plotKs(trPro50,ksPro50,obsPro50,1,26,false,50);
sgtitle('Prochlorococcus (Lagrangian, 10 dbar bin, Threshold = 50)');
exportgraphics(ax10a,'figures/ks_ProLag50.png'); clear ax10a;

%% 4.j.i. Synechococcus

ax11 = figure;
plotKs(trSyn,ksSyn,obsSyn,1,26,false);
sgtitle('Synechococcus (Lagrangian, 10 dbar bin)');
exportgraphics(ax11,'figures/ks_SynLag.png'); clear ax11;

%% 4.j.ii. Synechococcus (threshold = 50)

ax11a = figure;
plotKs(trSyn50,ksSyn50,obsSyn50,1,26,false,50);
sgtitle('Synechococcus (Lagrangian, 10 dbar bin, Threshold = 50)');
exportgraphics(ax11a,'figures/ks_SynLag50.png'); clear ax11a;

%% 4.k.i. Picoeukaryotes

ax12 = figure;
plotKs(trPic,ksPic,obsPic,1,26,false);
sgtitle('Picoeukaryotes (Lagrangian, 10 dbar bin)');
exportgraphics(ax12,'figures/ks_PicLag.png'); clear ax12;

%% 4.k.ii. Picoeukaryotes (threshold = 50)

ax12a = figure;
plotKs(trPic50,ksPic50,obsPic50,1,26,false,50);
sgtitle('Picoeukaryotes (Lagrangian, 10 dbar bin, Threshold = 50)');
exportgraphics(ax12a,'figures/ks_PicLag50.png'); clear ax12a;