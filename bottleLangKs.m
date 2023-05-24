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

idHet = num2str(importdata('data\hetBac_05-21.txt').data(:,1));
pHet = importdata('data\hetBac_05-21.txt').data(:,4);
het = importdata('data\hetBac_05-21.txt').data(:,5);

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
[trChla,ksChla,obsChla,skChla,kuChla,sdChla,c95chla,muChla] = ksOfLagrangian(idChla,pChla,dcm,chla,669);                        % chla (regular method)
b = 1:5890;
[trChla07,ksChla07,obsChla07,skChla07,kuChla07,sdChla07,c95chla07,muChla07] = ksOfLagrangian(idChla(b,:),pChla(b),dcm(b,:),chla(b),407);  % chla (regular method), 1988-2007
[trHplc,ksHplc,obsHplc,skHplc,kuHplc,sdHplc,c95hplc,muHplc] = ksOfLagrangian(idHplc,pHplc,dcm,hplc,317);                        % chla (HPLC method)
[trCmo,ksCmo,obsCmo,skCmo,kuCmo,sdCmo,c95cmo,muCmo] = ksOfLagrangian(idCmo,pCmo,dcm,cmo,271);                              % HPLC Chlorophyll a (Monovinyl)
[trCdi,ksCdi,obsCdi,skCdi,kuCdi,sdCdi,c95cdi,muCdi] = ksOfLagrangian(idCdi,pCdi,dcm,cdi,271);                              % HPLC Chlorophyll a (Divinyl)
[trCar,ksCar,obsCar,skCar,kuCar,sdCar,c95car,muCar] = ksOfLagrangian(idCar,pCar,dcm,car,332);                              % Particulate Carbon
[trNit,ksNit,obsNit,skNit,kuNit,sdNit,c95nit,muNit] = ksOfLagrangian(idNit,pNit,dcm,nit,333);                              % Particulate Nitrogen
[trPho,ksPho,obsPho,skPho,kuPho,sdPho,c95pho,muPho] = ksOfLagrangian(idPho,pPho,dcm,pho,331);                              % Particulate Phosphorus
% [trPho11,ksPho11,obsPho11,skPho11,kuPho11,sdPho11,c95pho11,muPho11] = ksOfLagrangian(idPho11,pPho11,dcm,pho11,90);                  % Particulate Phosphorus (11-21)
[trAtp,ksAtp,obsAtp,skAtp,kuAtp,sdAtp,c95atp,muAtp] = ksOfLagrangian(idAtp,pAtp,dcm,atp,318);                              % ATP
[trHet,ksHet,obsHet,skHet,kuHet,sdHet,c95het,muHet] = ksOfLagrangian(idHet,pHet,dcm,het,147);                              % Heterotrophic Bacteria
[trPro,ksPro,obsPro,skPro,kuPro,sdPro,c95pro,muPro] = ksOfLagrangian(idPro,pPro,dcm,pro,147);                              % Prochlorococcus
[trSyn,ksSyn,obsSyn,skSyn,kuSyn,sdSyn,c95syn,muSyn] = ksOfLagrangian(idSyn,pSyn,dcm,syn,147);                              % Synechococcus
[trPic,ksPic,obsPic,skPic,kuPic,sdPic,c95pic,muPic] = ksOfLagrangian(idPic,pPic,dcm,pic,147);                              % Picoeukaryotes

% Lower Threshold (=50)
[trPro50,ksPro50,obsPro50,skPro50,kuPro50,sdPro50,c95pro50,muPro50] = ksOfLagrangian(idPro,pPro,dcm,pro,147,50);                     % Prochlorococcus
[trSyn50,ksSyn50,obsSyn50,skSyn50,kuSyn50,sdSyn50,c95syn50,muSyn50] = ksOfLagrangian(idSyn,pSyn,dcm,syn,147,50);                     % Synechococcus
[trPic50,ksPic50,obsPic50,skPic50,kuPic50,sdPic50,c95pic50,muPic50] = ksOfLagrangian(idPic,pPic,dcm,pic,147,50);                     % Picoeukaryotes

% Alternative thresholds 'A'
% Based on adding ~20% new observed points
% Not including regular method chla.
[trHplcA,ksHplcA,obsHplcA,skHplcA,kuHplcA] = ksOfLagrangian(idHplc,pHplc,dcm,hplc,317,76);                  % chla (HPLC method)
[trCmoA,ksCmoA,obsCmoA,skCmoA,kuCmoA] = ksOfLagrangian(idCmo,pCmo,dcm,cmo,271,71);                        % HPLC Chlorophyll a (Monovinyl)
[trCdiA,ksCdiA,obsCdiA,skCdiA,kuCdiA] = ksOfLagrangian(idCdi,pCdi,dcm,cdi,271,68);                        % HPLC Chlorophyll a (Divinyl)
[trCarA,ksCarA,obsCarA,skCarA,kuCarA] = ksOfLagrangian(idCar,pCar,dcm,car,332,76);                        % Particulate Carbon
[trNitA,ksNitA,obsNitA,skNitA,kuNitA] = ksOfLagrangian(idNit,pNit,dcm,nit,333,76);                        % Particulate Nitrogen
[trPhoA,ksPhoA,obsPhoA,skPhoA,kuPhoA] = ksOfLagrangian(idPho,pPho,dcm,pho,331,63);                        % Particulate Phosphorus
[trPho11,ksPho11,obsPho11,skPho11,kuPho11] = ksOfLagrangian(idPho11,pPho11,dcm,pho11,50);                 % Particulate Phosphorus (11-21)
[trAtpA,ksAtpA,obsAtpA,skAtpA,kuAtpA] = ksOfLagrangian(idAtp,pAtp,dcm,atp,318,70);                        % ATP
[trHetA,ksHetA,obsHetA,skHetA,kuHetA] = ksOfLagrangian(idHet,pHet,dcm,het,147,46);                        % Heterotrophic Bacteria
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
plotKs(trHetA,ksHetA,obsHetA,skHetA,kuHetA,1,26,false,46);
sgtitle('Heterotrophic Bacteria (Lagrangian, 10 dbar, Threshold = 46)');
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

%% 5.1 Visualise mean

norm = 5; logn = 6;

xax = 1:1:20;
ax14 = figure;
subplot(2,1,1)
scatter(muChla(:,norm),trChla,'DisplayName','[chl a] (fluo.)','Marker','+');
hold on
scatter(muHplc(:,norm),trHplc,'DisplayName','[chl a] (HPLC)','Marker','+');
scatter(muAtp(:,norm),trAtp,'DisplayName','ATP','Marker','+');
scatter(muCmo(:,norm),trCmo,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
scatter(muCdi(:,norm),trCdi,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
scatter(muCar(:,norm),trCar,'DisplayName','Particulate Carbon','Marker','+');
scatter(muNit(:,norm),trNit,'DisplayName','Particulate Nitrogen','Marker','+');
scatter(muPho(:,norm),trPho,'DisplayName','Particulate Phosphorus');
scatter(muHet(:,norm),trHet,'DisplayName','Heterotrophic Bacteria ');
scatter(muPro(:,norm),trPro,'DisplayName','Prochlorococcus');
scatter(muSyn(:,norm),trSyn,'DisplayName','Synechococcus');
hold off
xlim([0.995 1.00]);
legend('Location','bestoutside'); set(gca,'YDir','reverse');
h = title('$\frac{\mu_{mle}}{\mu_{data}}$ normal','Interpreter','latex');
h.FontSize=30;

% tmpY = 5:10:200;
subplot(2,1,2)
scatter(muChla(:,logn),trChla,'DisplayName','[chl a] (fluo.)','Marker','+');
hold on
scatter(muHplc(:,logn),trHplc,'DisplayName','[chl a] (HPLC)','Marker','+');
scatter(muAtp(:,logn),trAtp,'DisplayName','ATP','Marker','+');
scatter(muCmo(:,logn),trCmo,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
scatter(muCdi(:,logn),trCdi,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
scatter(muCar(:,logn),trCar,'DisplayName','Particulate Carbon','Marker','+');
scatter(muNit(:,logn),trNit,'DisplayName','Particulate Nitrogen','Marker','+');
scatter(muPho(:,logn),trPho,'DisplayName','Particulate Phosphorus');
scatter(muHet(:,logn),trHet,'DisplayName','Heterotrophic Bacteria ');
scatter(muPro(:,logn),trPro,'DisplayName','Prochlorococcus');
scatter(muSyn(:,logn),trSyn,'DisplayName','Synechococcus');
hold off
legend('Location','bestoutside'); set(gca,'YDir','reverse');
h = title('$\frac{\mu_{mle}}{\mu_{data}}$ lognormal','Interpreter','latex');
h.FontSize=30;

sgtitle('Data vs. MLE: Lagrangian Mean');
exportgraphics(ax14,'figures/meanCompLag.png'); clear ax14;


%% 5.2 Visualise STD
% By definition, the standard deviation of the Lagrangian should be zero,
% so really these figures shouldn't show anything.

norm = 5; logn = 6;

colormap("parula");

xax = 1:1:20;
ax15 = figure;
subplot(2,1,1)
scatter(sdChla(:,norm),trChla,'DisplayName','[chl a] (fluo.)','Marker','+');
hold on
scatter(sdHplc(:,norm),trHplc,'DisplayName','[chl a] (HPLC)','Marker','+');
scatter(sdAtp(:,norm),trAtp,'DisplayName','ATP','Marker','+');
scatter(sdCmo(:,norm),trCmo,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
scatter(sdCdi(:,norm),trCdi,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
scatter(sdCar(:,norm),trCar,'DisplayName','Particulate Carbon','Marker','+');
scatter(sdNit(:,norm),trNit,'DisplayName','Particulate Nitrogen','Marker','+');
scatter(sdPho(:,norm),trPho,'DisplayName','Particulate Phosphorus');
scatter(sdHet(:,norm),trHet,'DisplayName','Heterotrophic Bacteria ');
scatter(sdPro(:,norm),trPro,'DisplayName','Prochlorococcus');
scatter(sdSyn(:,norm),trSyn,'DisplayName','Synechococcus');
hold off
xlim([0.995 1.00]);
legend('Location','bestoutside'); set(gca,'YDir','reverse');
h = title('$\frac{\sigma_{mle}}{\sigma_{data}}$ normal','Interpreter','latex');
h.FontSize=30;

% tmpY = 5:10:200;
subplot(2,1,2)
scatter(sdChla(:,logn),trChla,'DisplayName','[chl a] (fluo.)','Marker','+');
hold on
scatter(sdHplc(:,logn),trHplc,'DisplayName','[chl a] (HPLC)','Marker','+');
scatter(sdAtp(:,logn),trAtp,'DisplayName','ATP','Marker','+');
scatter(sdCmo(:,logn),trCmo,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
scatter(sdCdi(:,logn),trCdi,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
scatter(sdCar(:,logn),trCar,'DisplayName','Particulate Carbon','Marker','+');
scatter(sdNit(:,logn),trNit,'DisplayName','Particulate Nitrogen','Marker','+');
scatter(sdPho(:,logn),trPho,'DisplayName','Particulate Phosphorus');
scatter(sdHet(:,logn),trHet,'DisplayName','Heterotrophic Bacteria ');
scatter(sdPro(:,logn),trPro,'DisplayName','Prochlorococcus');
scatter(sdSyn(:,logn),trSyn,'DisplayName','Synechococcus');
hold off
xlim([0.995 1.00]);
legend('Location','bestoutside'); set(gca,'YDir','reverse');
h = title('$\frac{\sigma_{mle}}{\sigma_{data}}$ lognormal','Interpreter','latex');
h.FontSize=30;

sgtitle('Data vs. MLE: Lagrangian STD');
exportgraphics(ax15,'figures/stdCompLag.png'); clear ax15;


%% Visualise mean and error

meanMle = 1;

ax16 = figure;
subplot(2,2,1)
errorbar(muChla(:,meanMle),trChla,muChla(:,meanMle)-c95chla(:,1),c95chla(:,2)-muChla(:,meanMle),'horizontal');
set(gca,'YDir','reverse');
xlabel('ug/kg'); title('Fluoro Chl a');

subplot(2,2,2)
errorbar(muHplc(:,meanMle),trHplc,muHplc(:,meanMle)-c95hplc(:,1),c95hplc(:,2)-muHplc(:,meanMle),'horizontal');
set(gca,'YDir','reverse');
xlabel('ng/kg'); title('HPLC Chl a');

subplot(2,2,3)
errorbar(muAtp(:,meanMle),trAtp,muAtp(:,meanMle)-c95atp(:,1),c95atp(:,2)-muAtp(:,meanMle),'horizontal');
set(gca,'YDir','reverse');
xlabel('units'); title('ATP');

subplot(2,2,4)
errorbar(muCmo(:,meanMle),trCmo,muCmo(:,meanMle)-c95cmo(:,1),c95cmo(:,2)-muCmo(:,meanMle),'horizontal');
set(gca,'YDir','reverse');
xlabel('units'); title('Monovinyl Chl a');

exportgraphics(ax16,'figures/stdCompLag.png'); clear ax16;


% scatter(muChla(:,meanMle),trChla,'DisplayName','[chl a] (fluo.)','Marker','+');
% hold on
% scatter(c95chla(:,1),trChla);
% scatter(c95chla(:,2),trChla);

% scatter(muHplc(:,meanMle),trHplc,'DisplayName','[chl a] (HPLC)','Marker','+');
% scatter(muAtp(:,meanMle),trAtp,'DisplayName','ATP','Marker','+');
% scatter(muCmo(:,meanMle),trCmo,'DisplayName','[chl a] (HPLC Monovinyl)','Marker','+');
% scatter(muCdi(:,meanMle),trCdi,'DisplayName','[chl a] (HPLC Divinyl)','Marker','+');
% scatter(muCar(:,meanMle),trCar,'DisplayName','Particulate Carbon','Marker','+');
% scatter(muNit(:,meanMle),trNit,'DisplayName','Particulate Nitrogen','Marker','+');
% scatter(muPho(:,meanMle),trPho,'DisplayName','Particulate Phosphorus');
% scatter(muHet(:,meanMle),trHet,'DisplayName','Heterotrophic Bacteria ');
% scatter(muPro(:,meanMle),trPro,'DisplayName','Prochlorococcus');
% scatter(muSyn(:,meanMle),trSyn,'DisplayName','Synechococcus');
% hold off

% xlim([0.995 1.00]);
% legend('Location','best');
% % set(gca,'YDir','reverse');
% h = title('$\sigma_{mle}$ normal','Interpreter','latex');
% h.FontSize=15;
