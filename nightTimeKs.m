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

%%

% The following figures show the KS test results for single cruise-averaged
% fluorescence profiles across time
% I have NOT removed the individual suspicious casts from the averaging, so
% 

ax = figure;

% Isopycnal Eulerian
subplot(1,4,1)
plot(ksEi(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEi(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEi(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEi(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend('Location','best');
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Isopycnal Eulerian');

% Eulerian
subplot(1,4,2)
plot(ksEp(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEp(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEp(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEp(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian');

% Isopycnal Lagrangian
subplot(1,4,3)
plot(ksLi(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLi(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLi(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLi(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Isopycnal Lagrangian');

% Lagrangian
subplot(1,4,4)
plot(ksLp(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLp(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLp(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLp(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Lagrangian');

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

%% NIGHT
ax2 = figure;

% Isopycnal Eulerian
subplot(1,4,1)
plot(ksEiN(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEiN(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEiN(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEiN(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend('Location','best');
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Isopycnal Eulerian');

% Eulerian
subplot(1,4,2)
plot(ksEpN(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEpN(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpN(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEpN(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian');

% Isopycnal Lagrangian
subplot(1,4,3)
plot(ksLiN(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLiN(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLiN(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLiN(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Isopycnal Lagrangian');

% Lagrangian
subplot(1,4,4)
plot(ksLpN(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLpN(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLpN(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLpN(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Lagrangian');

sgtitle('Kolmogorov Smirnov: Night Time Cruise Average, 89-21');
exportgraphics(ax2,'figures/ks_newInclIsopycNight.png');
clear ax2;

%% NIGHT: 89-01 and 01-21 Comparison

% Cruise documentation states that HOT-130 was the last cruise with the
% original fluorometer; therefore we need to compare the ranges 1:130 and
% 131:329.

a = 1:130; b = 131:329;

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

clear a b;
%% Figures: Night 89-01
ax2a = figure;

% Isopycnal Eulerian
subplot(1,4,1)
plot(ksEiN_f1(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEiN_f1(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEiN_f1(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEiN_f1(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend('Location','best');
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Isopycnal Eulerian');

% Eulerian
subplot(1,4,2)
plot(ksEpN_f1(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEpN_f1(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpN_f1(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEpN_f1(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian');

% Isopycnal Lagrangian
subplot(1,4,3)
plot(ksLiN_f1(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLiN_f1(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLiN_f1(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLiN_f1(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Isopycnal Lagrangian');

% Lagrangian
subplot(1,4,4)
plot(ksLpN_f1(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLpN_f1(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLpN_f1(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLpN_f1(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Lagrangian');

sgtitle('Kolmogorov Smirnov: Night F1 89-01');
exportgraphics(ax2a,'figures/ks_yearly_f1.png');
clear ax2a;

%% Figures: Night 01-21
ax2b = figure;

% Isopycnal Eulerian
subplot(1,4,1)
plot(ksEiN_f2(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEiN_f2(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEiN_f2(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEiN_f2(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend('Location','best');
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Isopycnal Eulerian');

% Eulerian
subplot(1,4,2)
plot(ksEpN_f2(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEpN_f2(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpN_f2(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEpN_f2(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian');

% Isopycnal Lagrangian
subplot(1,4,3)
plot(ksLiN_f2(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLiN_f2(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLiN_f2(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLiN_f2(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Isopycnal Lagrangian');

% Lagrangian
subplot(1,4,4)
plot(ksLpN_f2(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLpN_f2(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLpN_f2(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLpN_f2(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Lagrangian');

sgtitle('Kolmogorov Smirnov: Night F2 01-21');
exportgraphics(ax2b,'figures/ks_yearly_f2.png');
clear ax2b;

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
% ksLiNa = getKS(meanIsoLFN,129,autVals);     % Isopycnal Lagrangian
ksLpNa = getKS(meanLpN,129,autVals);        % Lagrangian

%% SEASONAL KS FIGURES (night time)

%% Eulerian/Isopycnal/Night (EiN)

ax3 = figure;

% WINTER
subplot(1,4,1)
plot(ksEiNw(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEiNw(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEiNw(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEiNw(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Winter');

% SPRING
subplot(1,4,2)
plot(ksEiNsp(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEiNsp(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEiNsp(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEiNsp(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Spring');

% SUMMER
subplot(1,4,3)
plot(ksEiNsu(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEiNsu(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEiNsu(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEiNsu(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
h_leg = legend('Location','best');
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([0 250]);
ylabel('Pressure [db]');
title('Summer');

% AUTUMN
subplot(1,4,4)
plot(ksEiNa(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEiNa(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEiNa(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEiNa(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([0 250]);
ylabel('Pressure [db]');
title('Autumn');

% Make legend transparent
h_leg.BoxFace.ColorType='truecoloralpha';
h_leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');

sgtitle('Kolmogorov Smirnov: Eulerian Isopycnal Night, 89-21');
exportgraphics(ax3,'figures/ks_seasonal_EiN.png');
clear ax3;

%% Eulerian/Pressure/Night (EpN)

ax4 = figure;

% WINTER
subplot(1,4,1)
plot(ksEpNw(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEpNw(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpNw(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEpNw(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Winter');

% SPRING
subplot(1,4,2)
plot(ksEpNsp(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEpNsp(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpNsp(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEpNsp(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Spring');

% SUMMER
subplot(1,4,3)
plot(ksEpNsu(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEpNsu(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpNsu(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEpNsu(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
h_leg = legend('Location','best');
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([0 250]);
ylabel('Pressure [db]');
title('Summer');

% AUTUMN
subplot(1,4,4)
plot(ksEpNa(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEpNa(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEpNa(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEpNa(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([0 250]);
ylabel('Pressure [db]');
title('Autumn');

% Make legend transparent
h_leg.BoxFace.ColorType='truecoloralpha';
h_leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');

sgtitle('Kolmogorov Smirnov: Eulerian (Pressure Coord) Night, 89-21');
exportgraphics(ax4,'figures/ks_seasonal_EpN.png');
clear ax4;

%% Lagrangian/Isopycnal/Night (LiN)

ax5 = figure;

% WINTER
subplot(1,4,1)
plot(ksLiNw(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLiNw(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLiNw(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLiNw(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([-125 125]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Winter');

% SPRING
subplot(1,4,2)
plot(ksLiNsp(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLiNsp(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLiNsp(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLiNsp(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([-125 125]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Spring');

% SUMMER
subplot(1,4,3)
plot(ksLiNsu(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLiNsu(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLiNsu(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLiNsu(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
h_leg = legend('Location','best');
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Summer');

% AUTUMN - error calculating this, check back later
% subplot(1,4,4)
% plot(ksLiNa(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksEpNa(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksEpNa(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksEpNa(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylim([-125 125]);
% ylabel('Pressure [db]');
% title('Autumn');

% Make legend transparent
h_leg.BoxFace.ColorType='truecoloralpha';
h_leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');

sgtitle('Kolmogorov Smirnov: Lagrangian (Isopycnal Avg) Night, 89-21');
exportgraphics(ax5,'figures/ks_seasonal_LiN.png');
clear ax5;

%% Lagrangian/Pressure/Night (LpN)

ax6 = figure;

% WINTER
subplot(1,4,1)
plot(ksLpNw(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLpNw(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLpNw(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLpNw(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([-125 125]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Winter');

% SPRING
subplot(1,4,2)
plot(ksLpNsp(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLpNsp(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLpNsp(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLpNsp(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([-125 125]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Spring');

% SUMMER
subplot(1,4,3)
plot(ksLpNsu(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLpNsu(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLpNsu(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLpNsu(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
h_leg = legend('Location','best');
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Summer');

subplot(1,4,4)
plot(ksLpNa(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLpNa(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLpNa(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLpNa(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Autumn');

% Make legend transparent
h_leg.BoxFace.ColorType='truecoloralpha';
h_leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');

sgtitle('Kolmogorov Smirnov: Lagrangian (Pressure Avg) Night, 89-21');
exportgraphics(ax6,'figures/ks_seasonal_LpN.png');
clear ax6 h_leg;

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