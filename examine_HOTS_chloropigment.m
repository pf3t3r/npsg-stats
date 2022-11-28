% File to open and examine HOTS chloropigment data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',18);

%% Open data file and extract variables
data = importdata('data/hots-chloropigment.txt').data;

crn = data(:,1);
day = data(:,2);
pressure = data(:,3);
chloro = data(:,4);

chloro(chloro==-9) = NaN;

%% Examine First Days of Data

ax1 = figure;
plot(chloro(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
hold on
plot(chloro(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
plot(chloro(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
hold off
legend();
ylabel('pressure [db]');
xlabel('chloropigment (fluorescence) [ug/L]');
title('CTD: chloropigment vs depth (HOTS 88/89)');

exportgraphics(ax1,'figures/ctd-day-1.png');

%% Chlorophyll Depth- and Time-Series (1988-2022)

chloro2D = reshape(chloro,101,[]);
days = day + datetime(1988,09,30);
days = reshape(days,101,[]);
pres = reshape(pressure,101,[]);

[t_grid,p_grid] = meshgrid(datenum(days(1,:)),pres(:,1));

ax2 = figure;
[a,b] = contourf(t_grid,-p_grid,chloro2D,'LineColor','auto');
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
xlabel('Time');
ylabel('Depth [db]');
title('DCM Time Series: 1988 - 2021');

exportgraphics(ax2,'figures/fluorescence-1988-2021.png');