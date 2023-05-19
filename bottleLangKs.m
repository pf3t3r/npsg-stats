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

idPho11 = num2str(importdata('data\parP_11-21.txt').data(:,1));
pPho11 = importdata('data\parP_11-21.txt').data(:,4);
pho11 = importdata('data\parP_11-21.txt').data(:,5);

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
[trChla,ksChla,obsChla,skChla,kuChla] = ksOfLagrangian(idChla,pChla,dcm,chla,669);                        % chla (regular method)
b = 1:5890;
[trChla07,ksChla07,obsChla07,skChla07,kuChla07] = ksOfLagrangian(idChla(b,:),pChla(b),dcm(b,:),chla(b),407);  % chla (regular method), 1988-2007
[trHplc,ksHplc,obsHplc,skHplc,kuHplc] = ksOfLagrangian(idHplc,pHplc,dcm,hplc,317);                        % chla (HPLC method)
[trCmo,ksCmo,obsCmo,skCmo,kuCmo] = ksOfLagrangian(idCmo,pCmo,dcm,cmo,271);                              % HPLC Chlorophyll a (Monovinyl)
[trCdi,ksCdi,obsCdi,skCdi,kuCdi] = ksOfLagrangian(idCdi,pCdi,dcm,cdi,271);                              % HPLC Chlorophyll a (Divinyl)
[trCar,ksCar,obsCar,skCar,kuCar] = ksOfLagrangian(idCar,pCar,dcm,car,332);                              % Particulate Carbon
[trNit,ksNit,obsNit,skNit,kuNit] = ksOfLagrangian(idNit,pNit,dcm,nit,333);                              % Particulate Nitrogen
[trPho,ksPho,obsPho,skPho,kuPho] = ksOfLagrangian(idPho,pPho,dcm,pho,331);                              % Particulate Phosphorus
[trPho11,ksPho11,obsPho11,skPho11,kuPho11] = ksOfLagrangian(idPho11,pPho11,dcm,pho11,90);                  % Particulate Phosphorus (11-21)
[trAtp,ksAtp,obsAtp,skAtp,kuAtp] = ksOfLagrangian(idAtp,pAtp,dcm,atp,318);                              % ATP
[trHet,ksHet,obsHet,skHet,kuHet] = ksOfLagrangian(idHet,pHet,dcm,het,245);                              % Heterotrophic Bacteria
[trPro,ksPro,obsPro,skPro,kuPro] = ksOfLagrangian(idPro,pPro,dcm,pro,147);                              % Prochlorococcus
[trSyn,ksSyn,obsSyn,skSyn,kuSyn] = ksOfLagrangian(idSyn,pSyn,dcm,syn,147);                              % Synechococcus
[trPic,ksPic,obsPic,skPic,kuPic] = ksOfLagrangian(idPic,pPic,dcm,pic,147);                              % Picoeukaryotes

% Lower Threshold (=50)
[trPro50,ksPro50,obsPro50,skPro50,kuPro50] = ksOfLagrangian(idPro,pPro,dcm,pro,147,50);                     % Prochlorococcus
[trSyn50,ksSyn50,obsSyn50,skSyn50,kuSyn50] = ksOfLagrangian(idSyn,pSyn,dcm,syn,147,50);                     % Synechococcus
[trPic50,ksPic50,obsPic50,skPic50,kuPic50] = ksOfLagrangian(idPic,pPic,dcm,pic,147,50);                     % Picoeukaryotes

% Alternative thresholds 'A'
% Based on adding ~20% new observed points
% Not including regular method chla.
[trHplcA,ksHplcA,obsHplcA,skHplcA,kuHplcA] = ksOfLagrangian(idHplc,pHplc,dcm,hplc,317,76);                  % chla (HPLC method)
[trCmoA,ksCmoA,obsCmoA,skCmoA,kuCmoA] = ksOfLagrangian(idCmo,pCmo,dcm,cmo,271,71);                        % HPLC Chlorophyll a (Monovinyl)
[trCdiA,ksCdiA,obsCdiA,skCdiA,kuCdiA] = ksOfLagrangian(idCdi,pCdi,dcm,cdi,271,68);                        % HPLC Chlorophyll a (Divinyl)
[trCarA,ksCarA,obsCarA,skCarA,kuCarA] = ksOfLagrangian(idCar,pCar,dcm,car,332,76);                        % Particulate Carbon
[trNitA,ksNitA,obsNitA,skNitA,kuNitA] = ksOfLagrangian(idNit,pNit,dcm,nit,333,76);                        % Particulate Nitrogen
[trPhoA,ksPhoA,obsPhoA,skPhoA,kuPhoA] = ksOfLagrangian(idPho,pPho,dcm,pho,331,63);                        % Particulate Phosphorus
% [trPho11,ksPho11,obsPho11] = ksOfLagrangian(idPho11,pPho11,dcm,pho11,90);                 % Particulate Phosphorus (11-21)
[trAtpA,ksAtpA,obsAtpA,skAtpA,kuAtpA] = ksOfLagrangian(idAtp,pAtp,dcm,atp,318,70);                        % ATP
[trHetA,ksHetA,obsHetA,skHetA,kuHetA] = ksOfLagrangian(idHet,pHet,dcm,het,245,78);                        % Heterotrophic Bacteria
[trProA,ksProA,obsProA,skProA,kuProA] = ksOfLagrangian(idPro,pPro,dcm,pro,147,46);                        % Prochlorococcus
[trSynA,ksSynA,obsSynA,skSynA,kuSynA] = ksOfLagrangian(idSyn,pSyn,dcm,syn,147,48);                        % Synechococcus
[trPicA,ksPicA,obsPicA,skPicA,kuPicA] = ksOfLagrangian(idPic,pPic,dcm,pic,147,46);                        % Picoeukaryotes

%% 4. Plot the KS Statistics for each Bottle Concentration Time Series
% All of the following results are DCM-centred (Lagrangian), binned to 10
% dbar depths, and comprise the full HOT time series from 1988 up to 2021
% (2022- was not ready at time of download) unless otherwise specified.

%% 4.a.i. Chlorophyll a (Regular Method)

ax1 = figure;
plotKs(trChla,ksChla,obsChla,skChla,kuChla,10,35,false);
sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle chl-a [10 dbar bins] (88-21)');
exportgraphics(ax1,'figures/ks_bottleLagrangian10db.png'); clear ax1;

%% 4.a.ii. Chlorophyll a (Regular Method): 88-07

ax2 = figure;
plotKs(trChla07,ksChla07,obsChla07,skChla07,kuChla07,10,35,false);
sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle chl-a [10 dbar bins] (88-07)');
exportgraphics(ax2,'figures/ks_bottleLagrangian10db07.png'); clear ax2;

%% 4.b. HPLC Chlorophyll a

ax3 = figure;
plotKs(trHplc,ksHplc,obsHplc,skHplc,kuHplc,1,26,false);
sgtitle('HPLC Method: [chl a] (Lagrangian, 10 dbar bin)');
exportgraphics(ax3,'figures/ks_HplcLag.png'); clear ax3;

% Alternative Threshold
ax3a = figure;
plotKs(trHplcA,ksHplcA,obsHplcA,skHplcA,kuHplcA,1,26,false,76);
sgtitle('HPLC Method: [chl a] (Lagrangian, 10 dbar, Threshold = 76)');
exportgraphics(ax3a,'figures/ks_HplcLagA.png'); clear ax3a;

%% 4.c. HPLC Monovinyl Chlorophyll a

ax4 = figure;
plotKs(trCmo,ksCmo,obsCmo,skCmo,kuCmo,1,26,false);
sgtitle('HPLC Monovinyl Chlorophyll a (Lagrangian, 10 dbar bin)');
exportgraphics(ax4,'figures/ks_CmoLag.png'); clear ax4;

% Alternative Threshold
ax4a = figure;
plotKs(trCmoA,ksCmoA,obsCmoA,skCmoA,kuCmoA,1,26,false,71);
sgtitle('HPLC Monovinyl Chlorophyll a (Lagrangian, 10 dbar, Threshold = 71)');
exportgraphics(ax4a,'figures/ks_CmoLagA.png'); clear ax4a;

%% 4.d. HPLC Divinyl Chlorophyll a

ax5 = figure;
plotKs(trCdi,ksCdi,obsCdi,skCdi,kuCdi,1,26,false);
sgtitle('HPLC Divinyl Chlorophyll a (Lagrangian, 10 dbar bin)');
exportgraphics(ax5,'figures/ks_CdiLag.png'); clear ax5;

% Alternative Threshold
ax5a = figure;
plotKs(trCdiA,ksCdiA,obsCdiA,skCdiA,kuCdiA,1,26,false,68);
sgtitle('HPLC Divinyl Chlorophyll a (Lagrangian, 10 dbar, Threshold = 68)');
exportgraphics(ax5a,'figures/ks_CdiLagA.png'); clear ax5a;

%% 4.e. Particulate Carbon

ax6 = figure;
plotKs(trCar,ksCar,obsCar,skCar,kuCar,3,28,false);
sgtitle('Particulate Carbon (Lagrangian, 10 dbar bin)');
exportgraphics(ax6,'figures/ks_CarLag.png'); clear ax6;

% Alternative Threshold
ax6a = figure;
plotKs(trCarA,ksCarA,obsCarA,skCarA,kuCarA,3,28,false,76);
sgtitle('Particulate Carbon (Lagrangian, 10 dbar, Threshold = 76)');
exportgraphics(ax6a,'figures/ks_CarLagA.png'); clear ax6a;

%% 4.f. Particulate Nitrogen

ax7 = figure;
plotKs(trNit,ksNit,obsNit,skNit,kuNit,3,28,false);
sgtitle('Particulate Nitrogen (Lagrangian, 10 dbar bin)');
exportgraphics(ax7,'figures/ks_NitLag.png'); clear ax7;

ax7a = figure;
plotKs(trNitA,ksNitA,obsNitA,skNitA,kuNitA,3,28,false,76);
sgtitle('Particulate Nitrogen (Lagrangian, 10 dbar, Threshold = 76)');
exportgraphics(ax7a,'figures/ks_NitLag.png'); clear ax7a;

%% 4.g. Particulate Phosphorus

ax8 = figure;
plotKs(trPho,ksPho,obsPho,skPho,kuPho,3,28,false);
sgtitle('Particulate Phosphorus (Lagrangian, 10 dbar bin)');
exportgraphics(ax8,'figures/ks_PhoLag.png'); clear ax8;

ax8a = figure;
plotKs(trPho11,ksPho11,obsPho11,skPho11,kuPho11,3,28,false);
sgtitle('Particulate Phosphorus (Lagrangian, 10 dbar bin, 2011-2021)');
exportgraphics(ax8a,'figures/ks_PhoLag11.png'); clear ax8a;

ax8b = figure;
plotKs(trPhoA,ksPhoA,obsPhoA,skPhoA,kuPhoA,3,28,false,63);
sgtitle('Particulate Phosphorus (Lagrangian, 10 dbar, Threshold = 63)');
exportgraphics(ax8b,'figures/ks_PhoLagA.png'); clear ax8b;

%% 4.h. Heterotrophic Bacteria

ax9 = figure;
plotKs(trHet,ksHet,obsHet,skHet,kuHet,1,26,false);
sgtitle('Heterotrophic Bacteria (Lagrangian, 10 dbar bin)');
exportgraphics(ax9,'figures/ks_HetLag.png'); clear ax9;

ax9a = figure;
plotKs(trHetA,ksHetA,obsHetA,skHetA,kuHetA,1,26,false,78);
sgtitle('Heterotrophic Bacteria (Lagrangian, 10 dbar, Threshold = 78)');
exportgraphics(ax9a,'figures/ks_HetLagA.png'); clear ax9a;

%% 4.i. Prochlorococcus

ax10 = figure;
plotKs(trPro,ksPro,obsPro,skPro,kuPro,1,26,false);
sgtitle('Prochlorococcus (Lagrangian, 10 dbar bin)');
exportgraphics(ax10,'figures/ks_ProLag.png'); clear ax10;

ax10a = figure;
plotKs(trPro50,ksPro50,obsPro50,skPro50,kuPro50,1,26,false,50);
sgtitle('Prochlorococcus (Lagrangian, 10 dbar bin, Threshold = 50)');
exportgraphics(ax10a,'figures/ks_ProLag50.png'); clear ax10a;

ax10b = figure;
plotKs(trProA,ksProA,obsProA,skProA,kuProA,1,26,false,46);
sgtitle('Prochlorococcus (Lagrangian, 10 dbar, Threshold = 46)');
exportgraphics(ax10b,'figures/ks_ProLagA.png'); clear ax10b;

%% 4.j. Synechococcus

ax11 = figure;
plotKs(trSyn,ksSyn,obsSyn,skSyn,kuSyn,1,26,false);
sgtitle('Synechococcus (Lagrangian, 10 dbar)');
exportgraphics(ax11,'figures/ks_SynLag.png'); clear ax11;

ax11a = figure;
plotKs(trSyn50,ksSyn50,obsSyn50,skSyn50,kuSyn50,1,26,false,50);
sgtitle('Synechococcus (Lagrangian, 10 dbar, Threshold = 50)');
exportgraphics(ax11a,'figures/ks_SynLag50.png'); clear ax11a;

ax11b = figure;
plotKs(trSynA,ksSynA,obsSynA,skSynA,kuSynA,1,26,false,48);
sgtitle('Synechococcus (Lagrangian, 10 dbar, Threshold = 48)');
exportgraphics(ax11b,'figures/ks_SynLagA.png'); clear ax11b;

%% 4.k. Picoeukaryotes

ax12 = figure;
plotKs(trPic,ksPic,obsPic,skPic,kuPic,1,26,false);
sgtitle('Picoeukaryotes (Lagrangian, 10 dbar)');
exportgraphics(ax12,'figures/ks_PicLag.png'); clear ax12;

ax12a = figure;
plotKs(trPic50,ksPic50,obsPic50,skPic50,kuPic50,1,26,false,50);
sgtitle('Picoeukaryotes (Lagrangian, 10 dbar, Threshold = 50)');
exportgraphics(ax12a,'figures/ks_PicLag50.png'); clear ax12a;

ax12b = figure;
plotKs(trPicA,ksPicA,obsPicA,skPicA,kuPicA,1,26,false,46);
sgtitle('Picoeukaryotes (Lagrangian, 10 dbar, Threshold = 46)');
exportgraphics(ax12b,'figures/ks_PicLagA.png'); clear ax12b;


%% 4.l. ATP

ax13 = figure;
plotKs(trAtp,ksAtp,obsAtp,skAtp,kuAtp,2,27,false);
sgtitle('ATP (Lagrangian, 10 dbar)');
exportgraphics(ax13,'figures/ks_AtpLag.png'); clear ax13;

ax13a = figure;
plotKs(trAtpA,ksAtpA,obsAtpA,skAtpA,kuAtpA,2,27,false,70);
sgtitle('ATP (Lagrangian, 10 dbar, Threshold = 70)');
exportgraphics(ax13a,'figures/ks_AtpLagA.png'); clear ax13a;