% Fit distributions to a sample dataset.

clear;clc;close all;

% Import sample data.
chla = importdata("data\L0\chlaCtd250dbar.txt").data(:,4);
chla(chla==-9 | chla<=0) = nan;

% Fit distributions to the sample data.
N = fitdist(chla,"Normal");
L = fitdist(chla,"Lognormal");
G = fitdist(chla,"Gamma");
W = fitdist(chla,"Weibull");
% B = fitdist(chla,"Burr");

% Plot distributions
figure;
plot(N); hold on
plot(L);
plot(G);
plot(W); 
% plot(B);
hold off
legend("Normal","Lognormal","Gamma","Weibull");
