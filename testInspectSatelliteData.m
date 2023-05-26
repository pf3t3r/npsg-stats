%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
set(0,'defaultAxesFontSize',10);

%% Load Satellite Chl

gcaloha = load('data\aloha_chl_climatology.mat').gcaloha;

date = datetime(gcaloha.date,'ConvertFrom','datenum');
% date = gcaloha.date;
lon = gcaloha.lon;
lat = gcaloha.lat;
chl = gcaloha.chl;

%%

for i = 1:length(date)
    tmp = squeeze(chl(:,:,i));
    tmpMean(i) = mean(tmp,[1 2],'omitnan');
end

%%
figure
plot(date,tmpMean);
title('zonal and meridional mean chlorophyll from satellite data, station aloha');

%%

% Some extreme values before 2007, so let's try 2007-2021
% Alt: 2015 (6500) -> 2020 (8300)
% July 2011 (5040) - Sep 2012 (5470)
% up to July 2013 (5775)

tmp = tmpMean(5040:5775);
tmp(isnan(tmp)) = [];

tmp = movmean(tmp,3);

figure;
[mle,ks,~,~,~,~] = statsplot2(tmp);