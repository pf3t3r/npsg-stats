clear; clc; close all;
addpath("figures\kurtBias");
% This script will evaluate the performance of random data of varying
% sample size versus kurtosis and skewness. We aim to answer the following
% question: is there a bias for high kurtosis at low sample sizes?

%% Sample Size
s = 25000;

%% NORMAL case

kn = nan(100,1);
for i = 1:100
    a = randn(s,1);
    kn(i) = kurtosis(a);
end

%% LOGNORMAL case

kl = nan(100,1);
for i = 1:100
    b = lognrnd(1,0.3,[s 1]);
    kl(i) = kurtosis(b);
end

%% PLOT kn and kl
ax = figure;
plot(kn,DisplayName='Random Normal');
hold on
plot(kl,DisplayName='Random Lognormal');
hold off
legend();
title('Kurtosis of randomly-generated arrays',sprintf('Sample Size = %d',s));
str = "test" + s + "bias.png";
exportgraphics(ax,'figures/kurtBias/' + str); clear ax;