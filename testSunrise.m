close all; clc; clear;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 29 15]);
set(0,'defaultAxesFontSize',12);
addpath("baroneRoutines\");

%% Load data

f_e = load('datafiles\f_all.mat').f_er;
f_l = load('datafiles\f_all.mat').f_lag;
t = load('datafiles\f_all.mat').t_copy;

%% Find sunrise and sunset

lat = 22.75; lon = -158; date = t;
rs = nan(4347,2);

for i = 1:4347
    [tmp,~,~,~,~,~] = suncycle(lat,lon,date(i));
    tmp = tmp - 10;
    for j = 1:2
        if tmp(j) < 0
            tmp(j) = tmp(j) + 24;
        end
    end
    rs(i,:) = tmp;
end

%% Identify casts taken at night time

rs2 = hours(rs);
rs2.Format = 'hh:mm';

sunriseTime = hours(rs2(:,1)');
sunsetTime = hours(rs2(:,2)');
dayFrac = rem(datenum(t),1);
castTime = dayFrac*24;

castAtNight = nan(4347,1);

for i = 1:4347
    if castTime(i) < sunriseTime(i) || castTime(i) > sunsetTime(i)
        castAtNight(i) = 1;
    end
end

%% extract fE, fL, t for night time

nightCastIDs = [];

for i=1:4347
    if(~isnan(castAtNight(i)))
        disp(i);
        nightCastIDs = [nightCastIDs i];
    end
end

f_en = f_e(:,nightCastIDs);
f_ln = f_l(:,nightCastIDs);
t_n = t(nightCastIDs);


%% KS night time: calculate

% Eulerian
for i = 1:129
    disp(i);
    tmp = f_en(i,:);
    tmp(isnan(tmp)) = [];
    if length(tmp) > 2 
        [~,ks_en(:,i),~] = statsplot2(tmp,'noplot');
    end
end

% Lagrangian
for i = 1:129
    disp(i);
    tmp = f_ln(i,:);
    tmp(isnan(tmp)) = [];
    if length(tmp) > 2 
        [~,ks_ln(:,i),~] = statsplot2(tmp,'noplot');
    end
end

%% KS night time: figure

depthMeasurements = 129;
eulerianDepth = linspace(0,2*depthMeasurements,depthMeasurements);
lagrangianDepth = linspace(-128,128,depthMeasurements);

ax = figure;
ax.Position = [3 3 20 15];

% Eulerian
subplot(1,2,1)
plot(ks_en(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ks_en(2,:),eulerianDepth,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks_en(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks_en(4,:),eulerianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend();
ylim([0 250]);
set(gca,'YDir','reverse');
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian');

% Lagrangian
subplot(1,2,2)
plot(ks_ln(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ks_ln(2,:),lagrangianDepth,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks_ln(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks_ln(4,:),lagrangianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend();
ylim([-125 125]);
set(gca,'YDir','reverse');
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian');

sgtitle('Kolmogorov-Smirnov Test: night-time casts (1989-2021)');
exportgraphics(ax,'figures/ks_night_89-21.png');