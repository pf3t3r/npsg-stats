close all; clc; clear;

% File to open HOTS chloropigment data and examine it

data = importdata('data/hots-chloropigment.txt').data;

crn = data(:,1);
day = data(:,2);
pressure = data(:,3);
chloro = data(:,4);

ax = figure;
plot(chloro(1:101),-pressure(1:101),'DisplayName','First Day');
hold on
plot(chloro(102:202),-pressure(102:202),'DisplayName','Second Day');
hold off
legend();
ylabel('pressure [db]');
xlabel('chloropigment (fluorescence (?)) [ug/L]');
title('CTD: chloropigment vs depth (HOTS Day 1 - 1988)');

exportgraphics(ax,'figures/ctd-day-1.png');