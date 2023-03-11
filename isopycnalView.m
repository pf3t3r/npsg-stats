close all; clc; clear;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 29 15]);
set(0,'defaultAxesFontSize',11);
addpath("baroneRoutines\");

%%
iso = load('datafiles\ctd_iso_master2.mat').iso;
ctd = load('datafiles\ctd_iso_master2.mat').ctd;

cruises = 329;
missingCruises = [21,207];
cruisesRecorded = linspace(1,cruises,cruises);
cruisesRecorded(missingCruises) = [];
clear cruises missingCruises;

%%

sig = iso.sig;

f_iso = [];
t_iso = [];
pcm = [];
for i = cruisesRecorded
    tmpF = iso(i).f;
    tmpT = ctd(i).date;
    tmpPCM = ctd(i).pcm;
    f_iso = [f_iso tmpF];
    t_iso = [t_iso tmpT];
    pcm = [pcm tmpPCM];
end
clear tmpF tmpT i;
t_iso = datetime(t_iso,'ConvertFrom','datenum');

f_iso(f_iso < 0) = NaN;
depthMeasurements = 129;
eulerianDepth = linspace(0,2*depthMeasurements,depthMeasurements);
lagrangianDepth = linspace(-128,128,depthMeasurements);

%% Remove the following casts
% 1. All those casts which have NaN as the pressure at the DCM, and
% 2. All those casts that were flagged for having anomalously large mean
% values (>0.5 ug/kg) and/or large single values (>2 ug/kg), and were
% confirmed visually as suspicious.

% Step (1)
f_copy = f_iso;

p_copy = pcm;
p_copy(isnan(pcm)) = [];

f_copy(:,isnan(pcm)) = [];
f_isoR = f_copy(1:129,:);

t_copy = t_iso;
t_copy(isnan(pcm)) = [];

% Step (2)
dodgyCasts = [18 248 253 666 667 668 669 ...
    1213 1214 1215 1216 1217 1218 1219 1220 1221 1222 1223 1224 ...
    2124 2371 2434 2599 2610 2626 2898 2918 2926 2979 ...
    3110 3133 3560 3561];
f_isoR2 = f_isoR;
f_isoR2(:,dodgyCasts) = [];
t_copy(dodgyCasts) = [];

%% Find Lagrangian

[sig_cm,sig_cmID] = max(f_isoR2);
f_iso_lag = NaN(129,length(f_isoR2(1,:)));
offset_fISOL = 65 - sig_cmID;

for i = 1:length(f_isoR2(1,:))
    f_iso_lag(:,i) = circshift(f_isoR2(1:129,i),offset_fISOL(i));
    if offset_fISOL(i) > -1 && offset_fISOL(i) < 40
        disp(i);
        f_iso_lag(1:offset_fISOL(i),i) = NaN;
    elseif offset_fISOL(i) == -1
        f_iso_lag(end,i) = NaN;
    elseif offset_fISOL(i) < -1 && offset_fISOL(i) > -40
        disp(i);
        f_iso_lag((end+offset_fISOL(i)):end,i) = NaN;
    elseif abs(offset_fISOL(i)) > 40
        f_iso_lag(:,i) = NaN;
    end
end

%% ALL CASTS
%% Eul + Lag All Cast: Calculate KS

% Eulerian
for i = 1:129
    disp(i);
    tmp = f_isoR2(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksISO(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:129
    disp(i);
    tmp = f_iso_lag(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksISO_L(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;

%% Eul + Lag All Cast: KS Figs
sig_lag = sig(65) - sig(1:129);

ax = figure;
ax.Position = [3 3 13 15];

% Eulerian
subplot(1,2,1)
plot(ksISO(1,:),sig(1:129),'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksISO(2,:),sig(1:129),'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksISO(3,:),sig(1:129),'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksISO(4,:),sig(1:129),'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend();
xlim([0 0.8]);
set(gca,'YDir','reverse');
xlabel('p-value');
ylabel('\sigma^{\Theta} [kg m^{-3}]');
title('Eulerian');

% Lagrangian
subplot(1,2,2)
plot(ksISO_L(1,:),sig_lag,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksISO_L(2,:),sig_lag,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksISO_L(3,:),sig_lag,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksISO_L(4,:),sig_lag,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
% legend();
xlim([0 0.6]);
set(gca,'YDir','reverse');
xlabel('p-value');
ylabel('\Delta \sigma^{\Theta} [kg m^{-3}]');
title('Lagrangian');

sgtitle('KS on all casts (isopycnal, 89-21)');
exportgraphics(ax,'figures/ks_iso_allCast.png');

%% NIGHT CASTS only
%% Eul + Lag Night Cast: Calculate KS

nightCastID = load('datafiles\nightCast.mat').nightCastIDs;
f_iso_en = f_isoR2(:,nightCastID);
f_iso_ln = f_iso_lag(:,nightCastID);
% t_n = t(nightCastIDs);

% Eulerian
for i = 1:129
    disp(i);
    tmp = f_iso_en(i,:);
    tmp(isnan(tmp)) = [];
    if length(tmp) > 2
        [~,ksIsoEN(:,i),~] = statsplot2(tmp,'noplot');
    end
end
clear tmp;

% Lagrangian
for i = 1:129
    disp(i);
    tmp = f_iso_ln(i,:);
    tmp(isnan(tmp)) = [];
    if length(tmp) > 2
        [~,ksIsoLN(:,i),~] = statsplot2(tmp,'noplot');
    end
end
clear tmp;

%% Eul + Lag Night Casts: KS Figs

ax2 = figure;
ax2.Position = [3 3 13 15];

% Eulerian
subplot(1,2,1)
plot(ksIsoEN(1,:),sig(1:129),'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksIsoEN(2,:),sig(1:129),'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksIsoEN(3,:),sig(1:129),'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksIsoEN(4,:),sig(1:129),'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend();
xlim([0 0.8]);
set(gca,'YDir','reverse');
xlabel('p-value');
ylabel('\sigma^{\Theta} [kg m^{-3}]');
title('Eulerian');

% Lagrangian
subplot(1,2,2)
plot(ksIsoLN(1,:),sig_lag,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksIsoLN(2,:),sig_lag,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksIsoLN(3,:),sig_lag,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksIsoLN(4,:),sig_lag,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
% legend();
xlim([0 0.6]);
set(gca,'YDir','reverse');
xlabel('p-value');
ylabel('\Delta \sigma^{\Theta} [kg m^{-3}]');
title('Lagrangian');

sgtitle('KS on night casts (isopycnal, 89-21)');
exportgraphics(ax2,'figures/ks_iso_nightCast.png');

%% mean DCM

figure;
plot(mean(f_isoR2,2,'omitnan'),sig(1:129),'Color',[0.6 0.6 0.6]);
set(gca,'YDir','reverse');