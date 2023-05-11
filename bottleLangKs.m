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
end

% repeat crn
tcrn = repmat(crn',31,1)';

% Combine Above in One Array
tcast = reshape(cast',[],1);
ttcrn = reshape(tcrn',[],1);
tpcm = reshape(pcm',[],1);

dcmArray = [ttcrn tcast tpcm];
clear i tmp cast crn ctd pcm tcrn tcast ttcrn tpcm;

%% 2. Create Bottle Array

bottleId = num2str(importdata('data/hotbot-88_21.txt').data(:,1));
tCRN = str2num(bottleId(:,1:3)); %str2double doesn't work
tCast = str2num(bottleId(:,6:8));
tCast(tCast==100) = nan;
tP = importdata('data/hotbot-88_21.txt').data(:,4);

bottleArray = [tCRN tCast tP];
clear bottleId tCRN tCast tP;

%% Merge Arrays

tPcm = nan(9947,1);
tPLagrangian = nan(9947,1);

% find unique crn/cast combos and remove NaN cast indicators because it is
% unclear to which pcm they apply
t = rmmissing(unique(bottleArray(:,1:2),"rows"));

dcmArrayRowNo = []; x = 1;
for i = 1:10199
    if dcmArray(i,1:2) == t(x,1:2)
        dcmArrayRowNo = [dcmArrayRowNo i];
        x = x + 1;
    end
end

% save when bottleArray changes
tid = [];
for i = 2:9947
    if bottleArray(i,1) > bottleArray(i-1,1) || bottleArray(i,2) > bottleArray(i-1,2)
        tid = [tid i];
    end
end

tPcm(1:tid(1)-1) = dcmArray(dcmArrayRowNo(1),3);
for i = 2:669
    tPcm(tid(i):tid(i+1)-1) = dcmArray(dcmArrayRowNo(i),3);
end

bottleArray = [bottleArray tPcm];

tPLagrangian = bottleArray(:,3) - bottleArray(:,4);
bottleArray = [bottleArray tPLagrangian];

clear i x t dcmArrayRowNo tid tPcm tPLagrangian;

%% Bin pressure
pB10 = round(bottleArray(:,5),-1);
bottleArray = [bottleArray pB10];

%% Load chla
chla = importdata('data/hotbot-88_21.txt').data(:,5);
bottleArray = [bottleArray chla];

tmin = min(bottleArray(:,6));
tmax = max(bottleArray(:,6));
trange = tmin:10:tmax;

clear pB10 chla tmax tmin;
%% KS

chl = bottleArray(:,7);
pB = bottleArray(:,6);
ks = nan(5,41);
obsPerBin = nan(1,41);
for i = 1:length(trange)
    tmp = chl(pB==trange(i));
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obsPerBin(i) = length(tmp);
    if length(tmp) > 3
        disp(i);
        [~,ks(:,i),~] = statsplot2(tmp,'noplot');
    end
end
clear tmp;

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

for i = 1:41
    if obsPerBin(i) < 100
        ks(:,i) = nan;
    end
end

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
%% Chlorophyll a (Regular Method): Lagrangian 10 dbar: 88-07

% time = importdata('data/hotbot-88_21.txt').data(:,2);
% clear time;
chl07 = chl(1:5890);
pB07 = pB(1:5890);
ks07 = nan(5,41);
obsPerBin07 = nan(1,41);

for i = 1:length(trange)
    tmp1 = chl07(pB07==trange(i));
    tmp1(tmp1<=0) = nan;
    tmp1(isnan(tmp1)) = [];
    obsPerBin07(i) = length(tmp1);
    if length(tmp1) > 3
        disp(i);
        [~,ks07(:,i),~] = statsplot2(tmp1,'noplot');
    end
end
clear tmp1;

for i = 1:38
    if obsPerBin07(i) < 100
        ks07(:,i) = nan;
    end
end


ax2 = figure;
subplot(1,2,1)
barh(obsPerBin07);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([10 35]);
set(gca,"YTick",1:1:41,"YTickLabel",-240:10:130);
title('No. of Observations');

subplot(1,2,2)
plot(ks07(1,:),trange,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ks07(2,:),trange,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks07(3,:),trange,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks07(4,:),trange,'r.--','DisplayName','Gamma','MarkerSize',4);
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
