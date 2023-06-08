% This file will apply the Kolmogorov-Smirnov test to time- and
% depth-series of chl a from CTD fluorometry and bottle data. We will use
% only data from night-time casts (to remove the effect of
% non-photochemical quenching) and from after 2001 (specifically, from HOT
% cruise 130, when the current fluorometer was first used).

close all; clc; clear;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 29 15]);
set(0,'defaultAxesFontSize',11);
addpath("baroneRoutines\");

%% Load variables

iso = load('datafiles\ctd_iso_master2.mat').iso;
ctd = load('datafiles\ctd_iso_master2.mat').ctd;

cruises = 329;
missingCruises = [21,207];
cruisesRecorded = linspace(1,cruises,cruises);
cruisesRecorded(missingCruises) = [];
clear cruises missingCruises;

mlds = load('datafiles\MLD.mat').MLD;
mldt = load('datafiles\MLD.mat').MLDt;

%% Find (Eulerian) cruise-averaged fluorescence

meanEi = [];
meanEp = [];
for i = cruisesRecorded
    tmpI = mean([iso(i).f],2,"omitnan");
    meanEi = [meanEi tmpI];
    tmpE = mean([ctd(i).f],2,"omitnan");
    meanEp = [meanEp tmpE];
end
clear tmpI tmpE;

% Remove Negative Values
meanEi(meanEi<0) = nan;
meanEp(meanEp<0) = nan;

% This figure does not have casts removed: I need to incorporate the
% following sections
figure; plot(meanEi,[iso.sig]); set(gca,'YDir','reverse');
xlabel('Mean Fluorescence'); ylabel('Isopycnal Layer');

%% Find cruise-averaged fluorescence in Lagrangian coordinates

for i = cruisesRecorded
    % isopycnal
    tmp = iso(i).f(1:129,:);
    isoL(i).f = tmp;
    % pressure
    tmp2 = ctd(i).f(1:129,:);
    ctdL(i).f = tmp2;
end

for j = cruisesRecorded
    for i = 1:length(isoL(j).f(1,:)) % same length ctd and iso
        % iso
        tmp = isoL(j).f(:,i);
        [~,isoL(j).dcmID(i)] = max([isoL(j).f(:,i)]);
        isoL(j).offset(i) = 65 - isoL(j).dcmID(i);
        if isnan(isoL(j).offset(i))
            disp('i');
        else
            isoL(j).fL(:,i) = convertLagrangian(tmp,isoL(j).offset(i),329);
        end
        isoL(j).fLMean = mean(isoL(j).fL,2,'omitnan');

        % ctd
        tmp2 = ctdL(j).f(:,i);
        [~,ctdL(j).dcmID(i)] = max([ctdL(j).f(:,i)]);
        ctdL(j).offset(i) = 65 - ctdL(j).dcmID(i);
        if isnan(ctdL(j).offset(i))
            disp('i');
        else
            ctdL(j).fL(:,i) = convertLagrangian(tmp,ctdL(j).offset(i),329);
        end
        ctdL(j).fLMean = mean(ctdL(j).fL,2,'omitnan');
    end
end

clear tmp2;

%% Calculate KS

meanLi = [isoL.fLMean];
meanLi(meanLi<=0) = nan;
meanLp = [ctdL.fLMean];
meanLp(meanLp<0) = nan;

ksEi = getKS(meanEi,129);       % Isopycnal Eulerian
ksEp = getKS(meanEp,129);       % Eulerian
ksLi = getKS(meanLi,129);       % Isopycnal Lagrangian
ksLp = getKS(meanLp,129);       % Lagrangian

%% Bin the means

% edges = 0:10:250;
% p = ctd(1).p(1:129,:);
% test2 = meanEp(1:129,:);
% test = discretize([ctd(1).p(1:129,:)],edges);

%%

ax = figure;
plotKsCtd("yrly",ksEi,ksEp,ksLi,ksLp,[ctd(1).p(1:129)],mlds);
sgtitle('Kolmogorov Smirnov test: Cruise Average, 89-21');
exportgraphics(ax,'figures/ks_newInclIsopyc.png');
clear ax;

%% Night Time

lat = 22.75; lon = -158; 
rs = nan(329,31,2);
rsMean = nan(329,2);

for i = 1:329
    datetimeStruct(i).date = datetime(ctd(i).date,'ConvertFrom','datenum');
end

for i = 1:329
    for j = 1:length(datetimeStruct(i).date)
        [tmp,~,~,~,~,~] = suncycle(lat,lon,datetimeStruct(i).date(j));
        tmp = tmp - 10;
        for k = 1:2
            if tmp(k) < 0
                tmp(k) = tmp(k) + 24;
            end
        end
        rs(i,j,:) = tmp;
    end
    rsMean(i,1) = mean(rs(i,:,1),'omitnan');
    rsMean(i,2) = mean(rs(i,:,2),'omitnan');
end

rs2 = hours(rsMean);
rs2.Format = 'hh:mm';

clear k rs rsMean;
%% cast at night or not

% reconvert to decimal hour for comparison
sunrise = hours(rs2(:,1));
sunset = hours(rs2(:,2));

for i = 1:329
    dayFrac = rem(datenum([ctd(i).date]),1);
    datetimeStruct(i).castTime = dayFrac*24;
end
clear dayFrac rs2;

datetimeStruct(i).nightID = [];
for i = cruisesRecorded
    for j = 1:length(datetimeStruct(i).date)
        if datetimeStruct(i).castTime(j) < sunrise(i)
            tmp = j;
            datetimeStruct(i).nightID = [datetimeStruct(i).nightID tmp];
        end
    end
end

meanEiN = nan(501,329);
meanEpN = nan(501,329);
meanLiN = nan(129,329);
meanLpN = nan(129,329);

tmp = cruisesRecorded(cruisesRecorded~=48);
cruisesCTD = tmp(tmp~=218);
cruisesISO = cruisesCTD(cruisesCTD~=276);

for i = cruisesISO
    tmpEN = ctd(i).f(:,datetimeStruct(i).nightID);
    tmpEN = mean(tmpEN,2,"omitnan");
    meanEpN(:,i) = tmpEN;
    tmpMldEpN = ctd(i).mld003(datetimeStruct(i).nightID);
    meanMldEpN = mean(tmpMldEpN);
    tmpDcm = ctd(i).pcm(datetimeStruct(i).nightID);
    tmpMeanDcm(i) = mean(tmpDcm,'omitnan');
    tmpLN = ctdL(i).fL(:,datetimeStruct(i).nightID);
    tmpLN = mean(tmpLN,2,'omitnan');
    meanLpN(:,i) = tmpLN;
end

for i = cruisesISO
    tmpIN = iso(i).f(:,datetimeStruct(i).nightID);
    tmpIN = mean(tmpIN,2,'omitnan');
    meanEiN(:,i) = tmpIN;
    tmpILN = isoL(i).fL(:,datetimeStruct(i).nightID);
    tmpILN = mean(tmpILN,2,'omitnan');
    meanLiN(:,i) = tmpILN;
    
end

meanEiN(meanEiN<=0) = nan;
meanEpN(meanEpN<=0) = nan;
meanLiN(meanLiN<=0) = nan;
meanLpN(meanLpN<=0) = nan;

clear j tmpEN tmpILN tmpIN tmpLN tmp;
save('datafiles\lagrangianData.mat','isoL');
%% NIGHT

ksEiN = getKS(meanEiN,129);     % Isopycnal Eulerian
ksEpN = getKS(meanEpN,129);     % Eulerian
ksLiN = getKS(meanLiN,129);     % Isopycnal Lagrangian
ksLpN = getKS(meanLpN,129);     % Lagrangian

%% mean dcm night

meanDcmYrly = mean(tmpMeanDcm,'omitnan');

%% NIGHT: 89-01 and 01-21 Comparison

% Cruise documentation states that HOT-130 was the last cruise with the
% original fluorometer; therefore we need to compare the ranges 1:130 and
% 131:329.

a = 1:130; b = 131:329;

for i = a
    tmpMldEpN_f1 = ctd(i).mld003(datetimeStruct(i).nightID);
    meanMldEpN_f1 = mean(tmpMldEpN_f1);
end

for i = b
    tmpMldEpN_f2 = ctd(i).mld003(datetimeStruct(i).nightID);
    meanMldEpN_f2 = mean(tmpMldEpN_f2);
end

% mean PCM for F2, night time only
% for i = 1:length(b) % no WTF, this is wrong!
%     tmpMeanPcm(i) = mean([ctd(i).pcm(datetimeStruct(i).nightID)],'omitnan');
%     tmpSTD(i) = std([ctd(i).pcm(datetimeStruct(i).nightID)]);
% end
% ttMeanPcm = round(mean(tmpMeanPcm,'omitnan'));
% ttSTD = mean(tmpSTD,'omitnan');

% above method for finding mean and STD for EiN is incorrect, use below
% EpN: Mean, STD
[tmp2,pId_EpN] = max(meanEpN);
pId_EpN(pId_EpN==1) = NaN;
stdDcmEpN = 2*std(pId_EpN,'omitnan');
meanDcmEpn = 2*(mean(pId_EpN,'omitnan'));

% EiN: Mean, STD
[tmp,pId_EiN] = max(meanEiN);
pId_EiN(pId_EiN==1) = NaN;
stdDcmEiN = 2*std(pId_EiN,'omitnan');
meanDcmEin = 2*(mean(pId_EiN,'omitnan'));
clear tmp;

% Mean fluorescence: F1 vs F2 (EI)
mmEiN_f1 = mean(meanEiN(:,a),2,'omitnan');
mmEiN_f2 = mean(meanEiN(:,b),2,'omitnan');
ax2aa = figure;
subplot(2,2,1)
plot(mmEiN_f1,0:2:1000,'DisplayName','F1'); hold on;
plot(mmEiN_f2,0:2:1000,'DisplayName','F2'); hold off; 
xlabel('[chl] [{\mu}gL^{-1}]'); ylabel('Pressure [dbar]');
set(gca,'YDir','reverse'); title('Eulerian (Isopycnal)');

mmEpN_f1 = mean(meanEpN(:,a),2,'omitnan');
mmEpN_f2 = mean(meanEpN(:,b),2,'omitnan');
subplot(2,2,2)
plot(mmEpN_f1,0:2:1000,'DisplayName','F1'); hold on;
plot(mmEpN_f2,0:2:1000,'DisplayName','F2'); hold off;
legend('Location','best');
xlabel('[chl] [{\mu}gL^{-1}]'); ylabel('Pressure [dbar]');
set(gca,'YDir','reverse'); title('Eulerian (Pressure)');

mmLiN_f1 = mean(meanLiN(:,a),2,'omitnan');
mmLiN_f2 = mean(meanLiN(:,b),2,'omitnan');
subplot(2,2,3)
plot(mmLiN_f1,-128:2:128,'DisplayName','F1'); hold on;
plot(mmLiN_f2,-128:2:128,'DisplayName','F2'); hold off;
xlabel('[chl] [{\mu}gL^{-1}]'); ylabel('Pressure [dbar]');
set(gca,'YDir','reverse'); title('Lagrangian (Isopycnal)');

mmLpN_f1 = mean(meanLpN(:,a),2,'omitnan');
mmLpN_f2 = mean(meanLpN(:,b),2,'omitnan');
subplot(2,2,4)
plot(mmLpN_f1,-128:2:128,'DisplayName','F1'); hold on;
plot(mmLpN_f2,-128:2:128,'DisplayName','F2'); hold off;
xlabel('[chl] [{\mu}gL^{-1}]'); ylabel('Pressure [dbar]');
set(gca,'YDir','reverse'); title('Lagrangian (Pressure)');

sgtitle('F1 vs F2: Average Fluorescence');
exportgraphics(ax2aa,'figures/avgFluo_f1Vsf2.png');

%% NIGHT

ax2 = figure;
plotKsCtd("yrly",ksEiN,ksEpN,ksLiN,ksLpN,[ctd(1).p(1:129)],meanMldEpN,meanDcmEpn,stdDcmEpN,meanDcmEin,stdDcmEiN);
sgtitle('Kolmogorov Smirnov: Night Time Cruise Average, 89-21');
exportgraphics(ax2,'figures/ksYrlyNight.png');
clear ax2;

%% F1 vs F2: Find KS

% Night 89-01 (fluorometer f1)
ksEiN_f1 = getKS(meanEiN(:,a),129);     % Isopycnal Eulerian
ksEpN_f1 = getKS(meanEpN(:,a),129);     % Eulerian
ksLiN_f1 = getKS(meanLiN(:,a),129);     % Isopycnal Lagrangian
ksLpN_f1 = getKS(meanLpN(:,a),129);     % Lagrangian

% Night 01-21 (fluorometer f2)
ksEiN_f2 = getKS(meanEiN(:,b),129);     % Isopycnal Eulerian
ksEpN_f2 = getKS(meanEpN(:,b),129);     % Eulerian
ksLiN_f2 = getKS(meanLiN(:,b),129);     % Isopycnal Lagrangian
ksLpN_f2 = getKS(meanLpN(:,b),129);     % Lagrangian

%% Figures: Night 89-01


% EpN F1: Mean, STD
[tmp,pId_EpN_f1] = max(meanEpN(:,a));
pId_EpN_f1(pId_EpN_f1==1) = NaN;
stdDcmEpN_f1 = 2*std(pId_EpN_f1,'omitnan');
meanDcmEpn_f1 = 2*(mean(pId_EpN_f1,'omitnan'));

% EiN F1: Mean, STD
[tmp,pId_EiN_f1] = max(meanEiN(:,a));
pId_EiN_f1(pId_EiN_f1==1) = NaN;
stdDcmEiN_f1 = 2*std(pId_EiN_f1,'omitnan');
meanDcmEiN_f1 = 2*(mean(pId_EiN_f1,'omitnan'));


ax2a = figure;
plotKsCtd("yrly",ksEiN_f1,ksEpN_f1,ksLiN_f1,ksLpN_f1,[ctd(1).p(1:129)],meanMldEpN_f1,meanDcmEpn_f1,stdDcmEpN_f1,meanDcmEiN_f1,stdDcmEiN_f1);
sgtitle('Kolmogorov Smirnov: Night F1 89-01');
exportgraphics(ax2a,'figures/ks_yearly_f1.png');
clear ax2a;

%% Figures: Night 01-21

% EpN F1: Mean, STD
[tmp,pId_EpN_f2] = max(meanEpN(:,b));
pId_EpN_f2(pId_EpN_f2==1) = NaN;
stdDcmEpN_f2 = 2*std(pId_EpN_f2,'omitnan');
meanDcmEpn_f2 = 2*(mean(pId_EpN_f2,'omitnan'));

% EiN F1: Mean, STD
[tmp,pId_EiN_f2] = max(meanEiN(:,b));
pId_EiN_f2(pId_EiN_f2==1) = NaN;
stdDcmEiN_f2 = 2*std(pId_EiN_f2,'omitnan');
meanDcmEiN_f2 = 2*(mean(pId_EiN_f2,'omitnan'));

ax2b = figure;
plotKsCtd("yrly",ksEiN_f2,ksEpN_f2,ksLiN_f2,ksLpN_f2,[ctd(1).p(1:129)],meanMldEpN_f2,meanDcmEpn_f2,stdDcmEpN_f2,meanDcmEiN_f2,stdDcmEiN_f2);
sgtitle('Kolmogorov Smirnov: Night F2 01-21');
exportgraphics(ax2b,'figures/ksYrlyf2.png');
clear ax2b;

clear a b;
%% SEASONAL

for i = cruisesRecorded
    tmp = datetime([ctd(i).date],'ConvertFrom','datenum');
    [ctdL(i).win,ctdL(i).spr,ctdL(i).sum,ctdL(i).aut] = splitBySeasons(tmp,length(tmp));
end

% find means per season
winVals = []; sprVals = []; sumVals = []; autVals = [];
for i = cruisesRecorded
    if ~isnan([ctdL(i).win])
        winVals = [winVals i];
    end
    if ~isnan([ctdL(i).spr])
        sprVals = [sprVals i];
    end
    if ~isnan([ctdL(i).sum])
        sumVals = [sumVals i];
    end
    if ~isnan([ctdL(i).aut])
        autVals = [autVals i];
    end
end

clear i tmp;
save('datafiles\lagrangianData.mat','ctdL','-append');

%% KS NIGHT Seasonal
% These routines may take ~10-20s per season.

% MLD by season
meanMldWin = mean(mlds(winVals));
meanMldSpr = mean(mlds(sprVals));
meanMldSum = mean(mlds(sumVals));
meanMldAut = mean(mlds(autVals));

% DCM by season
meanDcmWin = mean(tmpMeanDcm(winVals),'omitnan');
meanDcmSpr = mean(tmpMeanDcm(sprVals),'omitnan');
meanDcmSum = mean(tmpMeanDcm(sumVals),'omitnan');
meanDcmAut = mean(tmpMeanDcm(autVals),'omitnan');

% DCM by season ("isopycnal")
meanDcmWinSig = round(2*mean(pId_EiN(winVals),'omitnan'));
meanDcmSprSig = round(2*mean(pId_EiN(sprVals),'omitnan'));
meanDcmSumSig = round(2*mean(pId_EiN(sumVals),'omitnan'));
meanDcmAutSig = round(2*mean(pId_EiN(autVals),'omitnan'));

% WINTER
ksEiNw = getKS(meanEiN,129,winVals);    % Isopycnal Eulerian
ksEpNw = getKS(meanEpN,129,winVals);    % Eulerian
ksLiNw = getKS(meanLiN,129,winVals);  % Isopycnal Lagrangian
ksLpNw = getKS(meanLpN,129,winVals);      % Lagrangian

% SPRING
ksEiNsp = getKS(meanEiN,129,sprVals);      % Isopycnal Eulerian
ksEpNsp = getKS(meanEpN,129,sprVals);      % Eulerian
ksLiNsp = getKS(meanLiN,129,sprVals);     % Isopycnal Lagrangian
ksLpNsp = getKS(meanLpN,129,sprVals);        % Lagrangian

% SUMMER
ksEiNsu = getKS(meanEiN,129,sumVals);      % Isopycnal Eulerian
ksEpNsu = getKS(meanEpN,129,sumVals);      % Eulerian
ksLiNsu = getKS(meanLiN,129,sumVals);     % Isopycnal Lagrangian
ksLpNsu = getKS(meanLpN,129,sumVals);        % Lagrangian

% AUTUMN
ksEiNa = getKS(meanEiN,129,autVals);      % Isopycnal Eulerian
ksEpNa = getKS(meanEpN,129,autVals);      % Eulerian
ksLiNa = getKS(meanLiN,129,autVals);     % Isopycnal Lagrangian
ksLpNa = getKS(meanLpN,129,autVals);        % Lagrangian

%% SEASONAL KS FIGURES (night time)


mld = [meanMldWin meanMldSpr meanMldSum meanMldAut];
dcm = [meanDcmWin meanDcmSpr meanDcmSum meanDcmAut];

dcmS = [meanDcmWinSig meanDcmSprSig meanDcmSumSig meanDcmAutSig];


%% Eulerian/Isopycnal/Night (EiN)

ax3 = figure;
plotKsCtd("seas",ksEiNw,ksEiNsp,ksEiNsu,ksEiNa,[ctd(1).p(1:129)],mld,dcmS);
sgtitle('Kolmogorov Smirnov: Eulerian Isopycnal Night, 89-21');
exportgraphics(ax3,'figures/ks_seasonal_EiN.png');
clear ax3;

%% Eulerian/Pressure/Night (EpN)

ax4 = figure;
plotKsCtd("seas",ksEpNw,ksEpNsp,ksEpNsu,ksEpNa,[ctd(1).p(1:129)],mld,dcm);
sgtitle('Kolmogorov Smirnov: Eulerian (Pressure Coord) Night, 89-21');
exportgraphics(ax4,'figures/ks_seasonal_EpN.png');
clear ax4;

%% Lagrangian/Isopycnal/Night (LiN)

ax5 = figure;
plotKsCtd("seas",ksLiNw,ksLiNsp,ksLiNsu,ksLiNa,[ctd(1).p(1:129)]-129);
sgtitle('Kolmogorov Smirnov: Lagrangian (Isopycnal Avg) Night, 89-21');
exportgraphics(ax5,'figures/ks_seasonal_LiN.png');
clear ax5;

%% Lagrangian/Pressure/Night (LpN)

ax6 = figure;
plotKsCtd("seas",ksLpNw,ksLpNsp,ksLpNsu,ksLpNa,[ctd(1).p(1:129)]-129);
sgtitle('Kolmogorov Smirnov: Lagrangian (Pressure Avg) Night, 89-21');
exportgraphics(ax6,'figures/ks_seasonal_LpN.png');
clear ax6;

%% F2 ssnl
% Run seasonal analysis using the 'new' fluorometer F2 only, i.e. from 2001
% to 2021.

b = 131:329;

f2winVals = winVals - 131; f2sprVals = sprVals - 131;
f2sumVals = sumVals - 131; f2autVals = autVals - 131;
f2winVals(f2winVals<0) = [];
f2sprVals(f2sprVals<0) = [];
f2sumVals(f2sumVals<0) = [];
f2autVals(f2autVals<=0) = [];


% MLD by season
meanMldWinF2 = mean(mlds(winVals(36:end)));
meanMldSprF2 = mean(mlds(sprVals(33:end)));
meanMldSumF2 = mean(mlds(sumVals(33:end)));
meanMldAutF2 = mean(mlds(autVals(30:end)));

% DCM by season
meanDcmWinF2 = mean(tmpMeanDcm(winVals(36:end)),'omitnan');
meanDcmSprF2 = mean(tmpMeanDcm(sprVals(33:end)),'omitnan');
meanDcmSumF2 = mean(tmpMeanDcm(sumVals(33:end)),'omitnan');
meanDcmAutF2 = mean(tmpMeanDcm(autVals(30:end)),'omitnan');

% DCM by season ("isopycnal")
meanDcmWinSigF2 = round(2*mean(pId_EiN(winVals(36:end)),'omitnan'));
meanDcmSprSigF2 = round(2*mean(pId_EiN(sprVals(33:end)),'omitnan'));
meanDcmSumSigF2 = round(2*mean(pId_EiN(sumVals(33:end)),'omitnan'));
meanDcmAutSigF2 = round(2*mean(pId_EiN(autVals(30:end)),'omitnan'));

% WINTER
ksEiNwF2 = getKS(meanEiN(:,b),129,f2winVals);    % Eulerian Isopycnal
ksEpNwF2 = getKS(meanEpN(:,b),129,f2winVals);    % Eulerian Pressure
ksLiNwF2 = getKS(meanLiN(:,b),129,f2winVals);    % Lagrangian Isopycnal
ksLpNwF2 = getKS(meanLpN(:,b),129,f2winVals);    % Lagrangian Pressure

% SPRING
ksEiNspF2 = getKS(meanEiN(:,b),129,f2sprVals);      % Eulerian Isopycnal
ksEpNspF2 = getKS(meanEpN(:,b),129,f2sprVals);      % Eulerian Pressure
ksLiNspF2 = getKS(meanLiN(:,b),129,f2sprVals);      % Lagrangian Isopycnal
ksLpNspF2 = getKS(meanLpN(:,b),129,f2sprVals);      % Lagrangian Pressure

% SUMMER
ksEiNsuF2 = getKS(meanEiN(:,b),129,f2sumVals);      % Eulerian Isopycnal
ksEpNsuF2 = getKS(meanEpN(:,b),129,f2sumVals);      % Eulerian Pressure
ksLiNsuF2 = getKS(meanLiN(:,b),129,f2sumVals);      % Lagrangian Isopycnal
ksLpNsuF2 = getKS(meanLpN(:,b),129,f2sumVals);      % Lagrangian Pressure

% AUTUMN
ksEiNaF2 = getKS(meanEiN(:,b),129,f2autVals);      % Eulerian Isopycnal
ksEpNaF2 = getKS(meanEpN(:,b),129,f2autVals);      % Eulerian Pressure
ksLiNaF2 = getKS(meanLiN(:,b),129,f2autVals);     % Lagrangian Isopycnal
ksLpNaF2 = getKS(meanLpN(:,b),129,f2autVals);        % Lagrangian Pressure

mldF2 = [meanMldWinF2 meanMldSprF2 meanMldSumF2 meanMldAutF2];
dcmF2 = [meanDcmWinF2 meanDcmSprF2 meanDcmSumF2 meanDcmAutF2];
dcmSF2 = [meanDcmWinSigF2 meanDcmSprSigF2 meanDcmSumSigF2 meanDcmAutSigF2];

%% F2 Eulerian/Isopycnal/Night (EiN)

ax7 = figure;
plotKsCtd("seas",ksEiNwF2,ksEiNspF2,ksEiNsuF2,ksEiNaF2,[ctd(1).p(1:129)],mldF2,dcmSF2);
sgtitle('Kolmogorov Smirnov: Eulerian Isopycnal Night, 01-21');
exportgraphics(ax7,'figures/ks_F2_seasonal_EiN.png');
clear ax7;

%% F2 Eulerian/Pressure/Night (EpN)

ax8 = figure;
plotKsCtd("seas",ksEpNwF2,ksEpNspF2,ksEpNsuF2,ksEpNaF2,[ctd(1).p(1:129)],mldF2,dcmF2);
sgtitle('Kolmogorov Smirnov: Eulerian (Pressure Coord) Night, 01-21');
exportgraphics(ax8,'figures/ks_F2_seasonal_EpN.png');
clear ax8;

%% F2 Lagrangian/Isopycnal/Night (LiN)

ax9 = figure;
plotKsCtd("seas",ksLiNwF2,ksLiNspF2,ksLiNsuF2,ksLiNaF2,[ctd(1).p(1:129)]-129);
sgtitle('Kolmogorov Smirnov: Lagrangian (Isopycnal Avg) Night, 01-21');
exportgraphics(ax9,'figures/ks_F2_seasonal_LiN.png');
clear ax9;

%% Lagrangian/Pressure/Night (LpN)

ax10 = figure;
plotKsCtd("seas",ksLpNwF2,ksLpNspF2,ksLpNsuF2,ksLpNaF2,[ctd(1).p(1:129)]-129);
sgtitle('Kolmogorov Smirnov: Lagrangian (Pressure Avg) Night, 01-21');
exportgraphics(ax10,'figures/ks_F2_seasonal_LpN.png');
clear ax10;

%% MAYBE I should use the 'dodgy cast' removal part here ...
% 1. All those casts which have NaN as the pressure at the DCM, and
% 2. All those casts that were flagged for having anomalously large mean
% values (>0.5 ug/kg) and/or large single values (>2 ug/kg), and were
% confirmed visually as suspicious.

% Step (1)
% f_copy = f_iso;
% 
% p_copy = pcm;
% p_copy(isnan(pcm)) = [];
% 
% f_copy(:,isnan(pcm)) = [];
% f_isoR = f_copy(1:129,:);
% 
% t_copy = t_iso;
% t_copy(isnan(pcm)) = [];
% 
% % Step (2)
% dodgyCasts = [18 248 253 666 667 668 669 ...
%     1213 1214 1215 1216 1217 1218 1219 1220 1221 1222 1223 1224 ...
%     2124 2371 2434 2599 2610 2626 2898 2918 2926 2979 ...
%     3110 3133 3560 3561];
% f_isoR2 = f_isoR;
% f_isoR2(:,dodgyCasts) = [];
% t_copy(dodgyCasts) = [];

%% SECTIONS hereafter maybe don't use anymore

%% Find Lagrangian
% 
% [sig_cm,sig_cmID] = max(f_isoR2);
% f_iso_lag = NaN(129,length(f_isoR2(1,:)));
% offset_fISOL = 65 - sig_cmID;
% 
% for i = 1:length(f_isoR2(1,:))
%     f_iso_lag(:,i) = circshift(f_isoR2(1:129,i),offset_fISOL(i));
%     if offset_fISOL(i) > -1 && offset_fISOL(i) < 40
%         disp(i);
%         f_iso_lag(1:offset_fISOL(i),i) = NaN;
%     elseif offset_fISOL(i) == -1
%         f_iso_lag(end,i) = NaN;
%     elseif offset_fISOL(i) < -1 && offset_fISOL(i) > -40
%         disp(i);
%         f_iso_lag((end+offset_fISOL(i)):end,i) = NaN;
%     elseif abs(offset_fISOL(i)) > 40
%         f_iso_lag(:,i) = NaN;
%     end
% end

%% ALL CASTS
%% Eul + Lag All Cast: Calculate KS
% 
% % Eulerian
% for i = 1:129
%     disp(i);
%     tmp = f_isoR2(i,:);
%     tmp(isnan(tmp)) = [];
%     [~,ksISO(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;
% 
% % Lagrangian
% for i = 1:129
%     disp(i);
%     tmp = f_iso_lag(i,:);
%     tmp(isnan(tmp)) = [];
%     [~,ksISO_L(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;

%% Eul + Lag All Cast: KS Figs

% incorrect
% sig_lag = sig(65) - sig(1:129);
% 
% axX = figure;
% axX.Position = [3 3 13 15];
% 
% % Eulerian
% subplot(1,2,1)
% plot(ksISO(1,:),sig(1:129),'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksISO(2,:),sig(1:129),'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksISO(3,:),sig(1:129),'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksISO(4,:),sig(1:129),'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% legend();
% xlim([0 0.8]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('\sigma^{\Theta} [kg m^{-3}]');
% title('Eulerian');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksISO_L(1,:),sig_lag,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksISO_L(2,:),sig_lag,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksISO_L(3,:),sig_lag,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksISO_L(4,:),sig_lag,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% % legend();
% xlim([0 0.6]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('\Delta \sigma^{\Theta} [kg m^{-3}]');
% title('Lagrangian');
% 
% sgtitle('KS on all casts (isopycnal, 89-21)');
% exportgraphics(axX,'figures/ks_iso_allCast.png');
% clear axX;

%% NIGHT CASTS only
% %% Eul + Lag Night Cast: Calculate KS
% 
% nightCastID = load('datafiles\nightCast.mat').nightCastIDs;
% f_iso_en = f_isoR2(:,nightCastID);
% f_iso_ln = f_iso_lag(:,nightCastID);
% % t_n = t(nightCastIDs);
% 
% % Eulerian
% for i = 1:129
%     disp(i);
%     tmp = f_iso_en(i,:);
%     tmp(isnan(tmp)) = [];
%     if length(tmp) > 2
%         [~,ksIsoEN(:,i),~] = statsplot2(tmp,'noplot');
%     end
% end
% clear tmp;
% 
% % Lagrangian
% for i = 1:129
%     disp(i);
%     tmp = f_iso_ln(i,:);
%     tmp(isnan(tmp)) = [];
%     if length(tmp) > 2
%         [~,ksIsoLN(:,i),~] = statsplot2(tmp,'noplot');
%     end
% end
% clear tmp;

%% Eul + Lag Night Casts: KS Figs

% ax2 = figure;
% ax2.Position = [3 3 13 15];
% 
% % Eulerian
% subplot(1,2,1)
% plot(ksIsoEN(1,:),sig(1:129),'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksIsoEN(2,:),sig(1:129),'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksIsoEN(3,:),sig(1:129),'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksIsoEN(4,:),sig(1:129),'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% legend();
% xlim([0 0.8]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('\sigma^{\Theta} [kg m^{-3}]');
% title('Eulerian');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksIsoLN(1,:),sig_lag,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksIsoLN(2,:),sig_lag,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksIsoLN(3,:),sig_lag,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksIsoLN(4,:),sig_lag,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% % legend();
% xlim([0 0.6]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('\Delta \sigma^{\Theta} [kg m^{-3}]');
% title('Lagrangian');
% 
% sgtitle('KS on night casts (isopycnal, 89-21)');
% exportgraphics(ax2,'figures/ks_iso_nightCast.png');

%% mean DCM
% 
% figure;
% plot(mean(f_isoR2,2,'omitnan'),sig(1:129),'Color',[0.6 0.6 0.6]);
% set(gca,'YDir','reverse');