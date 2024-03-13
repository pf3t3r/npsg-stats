clc;clear;close all;

d1 = load("data\HOT-341_LADCP\s2c1_2hz.mat");
d2 = load("data\HOT-341_LADCP\s2c2_2hz.mat");
d3 = load("data\HOT-341_LADCP\s2c3_2hz.mat");

c1 = d1.c; p1 = d1.p;
c2 = d2.c; p2 = d2.p;
c3 = d3.c; p3 = d3.p;

figure;
scatter(c1,p1,[],'.');
hold on
scatter(c2,p2,[],'.');
scatter(c3,p3,[],'.'); hold off
ylim([10 300]);
set(gca,"YDir","reverse");