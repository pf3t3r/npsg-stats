%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
set(0,'defaultAxesFontSize',10);

%% Load Satellite Chl a

gcaloha = load('data\aloha_chl_climatology.mat').gcaloha;

date = datetime(gcaloha.date,'ConvertFrom','datenum');
% date = gcaloha.date;
lon = gcaloha.lon;
lat = gcaloha.lat;
chl = gcaloha.chl;

%% Inspect time series at each grid cell

figure
tiledlayout(12,12);

for i = 1:12
    for j = 1:12
        nexttile
        plot(squeeze(chl(i,j,:)));
        set(gca,'YTick',[]); set(gca,'XTick',[]);
    end
end

sgtitle('Chl a Surface Concentration: 1997-2022, All Grid Cells');

%% Find 'Regional' Time-Mean Satellite Chl a
% Here I mean the mean chl a over the entire 12x12 region across the time
% period 1997 - 2022. So in this respect regional is a combined meridional
% and zonal mean.

for i = 1:length(date)
    tmp = squeeze(chl(:,:,i));
    zonMerMeanChla(i) = mean(tmp,[1 2],'omitnan');
end
clear i tmp;

%% Satellite Chl a: Regional Time-Mean
figure
plot(date,zonMerMeanChla);
title('zonal and meridional mean chlorophyll from satellite data, station aloha');

%% Find Time-Mean Satellite Chl a for each grid cell (2D Map)
testMean = mean(chl,3,'omitnan');

%% Time-Mean Satellite Chl a over 12x12 Grid
nb = 100;

figure; contourf(testMean,linspace(0.0683,0.0694,nb),'LineColor','auto');
xticklabels(lon); yticklabels(lat); set(gca,'YDir','reverse');
xlabel('Longitude'); ylabel('Latitude');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = '[chl a] (ug/kg)';
title('[chl a]: mean surface concentration (1997-2022)');
clear c nb;
%% Find Kolmogorov-Smirnov Statistic: Regional Time-Mean Chl a

% Some extreme values before 2007, so let's try 2007-2021
% Alt: 2015 (6500) -> 2020 (8300)
% July 2011 (5040) - Sep 2012 (5470), or up to July 2013 (5775)

limA = 1;
limB = 9109;

randSam = datasample(squeeze(chl(1,1,:)),1000);
% tmp = zonMerMeanChla(limA:limB);
tmp = randSam;
tmp(isnan(tmp)) = [];

tmp = movmean(tmp,10);

figure;
[mle,ks,~,~,~,~] = statsplot2(tmp);

% clear tmp;