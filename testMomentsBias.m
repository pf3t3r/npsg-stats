clear; clc; close all;
addpath("figures\kurtBias");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

% This script will evaluate the performance of random data of varying
% sample size versus kurtosis and skewness. We aim to answer the following
% question: is there a bias for high kurtosis at low sample sizes? And is
% there a similar bias for skewness?

%% Kurtosis Bias

[kn1,kl1] = testKurtosisBiasHelper(16);
[kn2,kl2] = testKurtosisBiasHelper(32);
[kn3,kl3] = testKurtosisBiasHelper(64);
[kn4,kl4] = testKurtosisBiasHelper(128);
[kn5,kl5] = testKurtosisBiasHelper(256);
[kn6,kl6] = testKurtosisBiasHelper(512);
[kn7,kl7] = testKurtosisBiasHelper(1024);
[kn8,kl8] = testKurtosisBiasHelper(2048);
[kn9,kl9] = testKurtosisBiasHelper(4096);
[knX,klX] = testKurtosisBiasHelper(8192);
[kn11,kl11] = testKurtosisBiasHelper(16384);
[kn12,kl12] = testKurtosisBiasHelper(32768);
[kn13,kl13] = testKurtosisBiasHelper(65536);
[kn14,kl14] = testKurtosisBiasHelper(131072);

n1 = mean(kn1); l1 = mean(kl1);
n2 = mean(kn2); l2 = mean(kl2);
n3 = mean(kn3); l3 = mean(kl3);
n4 = mean(kn4); l4 = mean(kl4);
n5 = mean(kn5); l5 = mean(kl5);
n6 = mean(kn6); l6 = mean(kl6);
n7 = mean(kn7); l7 = mean(kl7);
n8 = mean(kn8); l8 = mean(kl8);
n9 = mean(kn9); l9 = mean(kl9);
nX = mean(knX); lX = mean(klX);
n11 = mean(kn11); l11 = mean(kl11);
n12 = mean(kn12); l12 = mean(kl12);
n13 = mean(kn13); l13 = mean(kl13);
n14 = mean(kn14); l14 = mean(kl14);

ax = figure;
subplot(1,2,1)
scatter(16,n1);
hold on
scatter(32,n2);
scatter(64,n3);
scatter(128,n4);
scatter(256,n5);
scatter(512,n6);
scatter(1028,n7);
scatter(2056,n8);
scatter(4128,n9);
scatter(8256,nX);
scatter(16384,n11);
scatter(32768,n12);
scatter(65536,n13);
scatter(131072,n14);
hold off
grid on;
ylim([2.1 3.9]);
set(gca,'xscale','log');
title('Normal');

subplot(1,2,2)
scatter(16,l1);
hold on
scatter(32,l2);
scatter(64,l3);
scatter(128,l4);
scatter(256,l5);
scatter(512,l6);
scatter(1028,l7);
scatter(2056,l8);
scatter(4128,l9);
scatter(8256,lX);
scatter(16384,l11);
scatter(32768,l12);
scatter(65536,l13);
scatter(131072,l14);
hold off
grid on;
set(gca,'xscale','log');
title('Lognormal');
sgtitle('Kurtosis Bias');

exportgraphics(ax,'figures/kurtBias/compOfBias.png');

clear;

%% Skewness Bias

[sn1,sl1] = testSkewnessBiasHelper(16);
[sn2,sl2] = testSkewnessBiasHelper(32);
[sn3,sl3] = testSkewnessBiasHelper(64);
[sn4,sl4] = testSkewnessBiasHelper(128);
[sn5,sl5] = testSkewnessBiasHelper(256);
[sn6,sl6] = testSkewnessBiasHelper(512);
[sn7,sl7] = testSkewnessBiasHelper(1024);
[sn8,sl8] = testSkewnessBiasHelper(2048);
[sn9,sl9] = testSkewnessBiasHelper(4096);
[snX,slX] = testSkewnessBiasHelper(8192);
[sn11,sl11] = testSkewnessBiasHelper(16384);
[sn12,sl12] = testSkewnessBiasHelper(32768);
[sn13,sl13] = testSkewnessBiasHelper(65536);
[sn14,sl14] = testSkewnessBiasHelper(131072);

n1 = mean(sn1); l1 = mean(sl1);
n2 = mean(sn2); l2 = mean(sl2);
n3 = mean(sn3); l3 = mean(sl3);
n4 = mean(sn4); l4 = mean(sl4);
n5 = mean(sn5); l5 = mean(sl5);
n6 = mean(sn6); l6 = mean(sl6);
n7 = mean(sn7); l7 = mean(sl7);
n8 = mean(sn8); l8 = mean(sl8);
n9 = mean(sn9); l9 = mean(sl9);
nX = mean(snX); lX = mean(slX);
n11 = mean(sn11); l11 = mean(sl11);
n12 = mean(sn12); l12 = mean(sl12);
n13 = mean(sn13); l13 = mean(sl13);
n14 = mean(sn14); l14 = mean(sl14);

ax = figure;
subplot(1,2,1)
scatter(16,n1);
hold on
scatter(32,n2);
scatter(64,n3);
scatter(128,n4);
scatter(256,n5);
scatter(512,n6);
scatter(1028,n7);
scatter(2056,n8);
scatter(4128,n9);
scatter(8256,nX);
scatter(16384,n11);
scatter(32768,n12);
scatter(65536,n13);
scatter(131072,n14);
hold off
grid on;
ylim([-0.2 0.2]);
set(gca,'xscale','log');
title('Normal');

subplot(1,2,2)
scatter(16,l1);
hold on
scatter(32,l2);
scatter(64,l3);
scatter(128,l4);
scatter(256,l5);
scatter(512,l6);
scatter(1028,l7);
scatter(2056,l8);
scatter(4128,l9);
scatter(8256,lX);
scatter(16384,l11);
scatter(32768,l12);
scatter(65536,l13);
scatter(131072,l14);
hold off
grid on;
set(gca,'xscale','log');
title('Lognormal');
sgtitle('Skewness Bias');

exportgraphics(ax,'figures/skewBias/compOfBias.png');