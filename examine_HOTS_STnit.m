% ALOHA (A Long-Term Oligotrophic Habitat Assessment; 22° 45'N, 158° 00'W)
% 22.75 N, 158 W

% File to open and examine HOTS T, S, Nitrate Data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',16);

%% Open data file and extract variables
data = importdata('data/hots-T-S-nit.txt').data;

crn = data(:,1);
day = data(:,2);
pressure = data(:,3);
T = data(:,4);
S = data(:,5);
nitrate = data(:,6);

T(T==-9) = NaN;
S(S==-9) = NaN;
nitrate(nitrate==-9) = NaN;

%% Calculate Conservative Temperature, Absolute Salinity, and Density

SA = gsw_SA_from_SP(S,pressure,-158,22.75);
CT = gsw_CT_from_t(SA,T,pressure);
rho = gsw_rho(SA,CT,pressure);

%% Examine First Days of Data for CT, SA, and sigma

ax1 = figure;
sgtitle('CTD: Temperature & Salinity vs Depth (HOTS 88/89)');
subplot(1,3,1)
plot(CT(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
hold on
plot(CT(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
plot(CT(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
legend('Location','best');
ylabel('pressure [db]');
xlabel('Conservative Temperature \Theta (Celsius)');

subplot(1,3,2)
plot(SA(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
hold on
plot(SA(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
plot(SA(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
hold off
legend('Location','best');
ylabel('pressure [db]');
xlabel('Absolute Salinity S_A (g/kg)');

subplot(1,3,3)
plot(rho(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
hold on
plot(rho(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
plot(rho(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
hold off
legend('Location','best');
ylabel('pressure [db]');
xlabel('Potential Density \sigma (kgm^{-3})');

exportgraphics(ax1,'figures/ctd-day-1-SA-CT-sigma.png');

%% SA-CT Plot

p_ref = 1000;
isopycs = 24:0.5:33;

% Adjust size because GSW preset gives a square SA-CT plot
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [16 5 20 20]);

ax2 = figure;
gsw_SA_CT_plot(SA,CT,p_ref,isopycs,'\it{S}\rm_A - {\Theta} diagram');
exportgraphics(ax2,'figures/ctd-SA-CT.png');