%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");
addpath("output\");
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

ax1 = figure;
tiledlayout(12,12);

for i = 1:12
    for j = 1:12
        nexttile
        plot(squeeze(chl(i,j,:)));
        set(gca,'YTick',[]); set(gca,'XTick',[]);
    end
end

sgtitle('Chl a Surface Concentration: 1997-2022, All Grid Cells');
exportgraphics(ax1,'figures/satellite/allCellTimeseries.png');

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
ax2 = figure;
plot(date,zonMerMeanChla);
title('zonal and meridional mean chlorophyll from satellite data, station aloha');
exportgraphics(ax2,'figures/satellite/zonMerMeanSurfaceChl.png');

%% Find Time-Mean Satellite Chl a for each grid cell (2D Map)
testMean = mean(chl,3,'omitnan');

%% Time-Mean Satellite Chl a over 12x12 Grid
nb = 100;

ax3 = figure; contourf(testMean,linspace(0.0683,0.0694,nb),'LineColor','auto');
xticklabels(lon); yticklabels(lat); set(gca,'YDir','reverse');
xlabel('Longitude'); ylabel('Latitude');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = '[chl a] (ug/kg)';
title('[chl a]: mean surface concentration (1997-2022)');
exportgraphics(ax3,'figures/satellite/chlContour.png');
clear c nb;
%% Find Kolmogorov-Smirnov Statistic: Regional Time-Mean Chl a

% Some extreme values before 2007, so let's try 2007-2021
% Alt: 2015 (6500) -> 2020 (8300)
% July 2011 (5040) - Sep 2012 (5470), or up to July 2013 (5775)

limA = 1;
limB = 9109;

randSam = datasample(squeeze(chl(1,1,:)),1000);
% tmp = zonMerMeanChla(limA:limB);
tmp = squeeze(chl(6,6,:));
% tmp = randSam;
tmp(isnan(tmp)) = [];

ax4 = figure;
[mle,ks,~] = statsplot2(tmp);
exportgraphics(ax4,'figures/satellite/hists.png');

%% Sensitivity Analysis

sampleNo = [100 250 500 1000 1500 2000];
movMean = [1 3 4 5 6 7 8];
mnKsS = nan(6,7,5);

if(exist("datafiles\pvalueSensitivityTest.mat","file"))
    disp('Loading...')
    load('pvalueSensitivityTest.mat', 'mnKsS');
else
    for i = 1:length(sampleNo)
        disp(i);
        for j = 1:length(movMean)
            for k = 1:10
                randSam = datasample(squeeze(chl(6,6,:)),sampleNo(i));
                randSam(isnan(randSam)) = [];
                randSam = movmean(randSam,movMean(j));
                [~,ksS(k,:),~,~,~,~] = statsplot2(randSam);
            end    
            mnKsS(i,j,:) = mean(ksS);
        end
    end
    save datafiles/pvalueSensitivityTest.mat mnKsS;
end


%% Sensitivity Analysis: Sample No. vs p-value (no moving mean)

ax5 = figure;
plot(sampleNo,mnKsS(:,1,1),'o-','Color','#a6cee3','DisplayName','normal');
hold on
plot(sampleNo,mnKsS(:,1,2),'+--','Color','#1f78b4','DisplayName','lognormal');
plot(sampleNo,mnKsS(:,1,3),'x-','Color','#b2df8a','DisplayName','weibull');
plot(sampleNo,mnKsS(:,1,4),'.--','Color','#33a02c','DisplayName','gamma');
hold off
legend();
ylabel('p-value'); xlabel('sample no.');
title('Station Aloha Surface Chl a: KS Test p-value vs. sample no.');
exportgraphics(ax5,'figures/satellite/ksVsP_sample.png');

%% Sensitivity Analysis: Including smoothing

ax6 = figure;
subplot(1,7,1)
plot(sampleNo,mnKsS(:,1,1),'o-','Color','#a6cee3','DisplayName','normal');
hold on
plot(sampleNo,mnKsS(:,1,2),'+--','Color','#1f78b4','DisplayName','lognormal');
plot(sampleNo,mnKsS(:,1,3),'x-','Color','#b2df8a','DisplayName','weibull');
plot(sampleNo,mnKsS(:,1,4),'.--','Color','#33a02c','DisplayName','gamma');
hold off
% legend();
ylabel('p-value'); xlabel('sample no.');
title('no smoothing');

subplot(1,7,2)
plot(sampleNo,mnKsS(:,2,1),'o-','Color','#a6cee3','DisplayName','normal');
hold on
plot(sampleNo,mnKsS(:,2,2),'+--','Color','#1f78b4','DisplayName','lognormal');
plot(sampleNo,mnKsS(:,2,3),'x-','Color','#b2df8a','DisplayName','weibull');
plot(sampleNo,mnKsS(:,2,4),'.--','Color','#33a02c','DisplayName','gamma');
hold off
% legend();
ylabel('p-value'); xlabel('sample no.');
title('smoothed','window = 3');

subplot(1,7,3)
plot(sampleNo,mnKsS(:,3,1),'o-','Color','#a6cee3','DisplayName','normal');
hold on
plot(sampleNo,mnKsS(:,3,2),'+--','Color','#1f78b4','DisplayName','lognormal');
plot(sampleNo,mnKsS(:,3,3),'x-','Color','#b2df8a','DisplayName','weibull');
plot(sampleNo,mnKsS(:,3,4),'.--','Color','#33a02c','DisplayName','gamma');
hold off
% legend();
ylabel('p-value'); xlabel('sample no.');
title('smoothed','window = 4')

subplot(1,7,4)
plot(sampleNo,mnKsS(:,4,1),'o-','Color','#a6cee3','DisplayName','normal');
hold on
plot(sampleNo,mnKsS(:,4,2),'+--','Color','#1f78b4','DisplayName','lognormal');
plot(sampleNo,mnKsS(:,4,3),'x-','Color','#b2df8a','DisplayName','weibull');
plot(sampleNo,mnKsS(:,4,4),'.--','Color','#33a02c','DisplayName','gamma');
hold off
% legend();
ylabel('p-value'); xlabel('sample no.');
title('smoothed','window = 5');

subplot(1,7,5)
plot(sampleNo,mnKsS(:,5,1),'o-','Color','#a6cee3','DisplayName','normal');
hold on
plot(sampleNo,mnKsS(:,5,2),'+--','Color','#1f78b4','DisplayName','lognormal');
plot(sampleNo,mnKsS(:,5,3),'x-','Color','#b2df8a','DisplayName','weibull');
plot(sampleNo,mnKsS(:,5,4),'.--','Color','#33a02c','DisplayName','gamma');
hold off
% legend();
ylabel('p-value'); xlabel('sample no.');
title('smoothed','window = 6');

subplot(1,7,6)
plot(sampleNo,mnKsS(:,6,1),'o-','Color','#a6cee3','DisplayName','normal');
hold on
plot(sampleNo,mnKsS(:,6,2),'+--','Color','#1f78b4','DisplayName','lognormal');
plot(sampleNo,mnKsS(:,6,3),'x-','Color','#b2df8a','DisplayName','weibull');
plot(sampleNo,mnKsS(:,6,4),'.--','Color','#33a02c','DisplayName','gamma');
hold off
% legend();
ylabel('p-value'); xlabel('sample no.');
title('smoothed','window = 7');

subplot(1,7,7)
plot(sampleNo,mnKsS(:,7,1),'o-','Color','#a6cee3','DisplayName','normal');
hold on
plot(sampleNo,mnKsS(:,7,2),'+--','Color','#1f78b4','DisplayName','lognormal');
plot(sampleNo,mnKsS(:,7,3),'x-','Color','#b2df8a','DisplayName','weibull');
plot(sampleNo,mnKsS(:,7,4),'.--','Color','#33a02c','DisplayName','gamma');
hold off
legend(); ylabel('p-value'); xlabel('sample no.');
title('smoothed','window = 8');

sgtitle('mean p-value vs. sample no.');
exportgraphics(ax6,'figures/satellite/ksVsPsmoothing.png');

% randSam = datasample(squeeze(chl(6,6,:)),1000);
