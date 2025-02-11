% Check relationship between sample size and skewness-kurtosis.
% This script will evaluate the performance of random data of varying
% sample size versus kurtosis and skewness. We aim to answer the following
% question: is there a bias for high kurtosis at low sample sizes? And is
% there a similar bias for skewness?
% Specifically, this script
% - creates random data sets of sample size s for four possible
% distributions [1]
% - r these multiple times and saves kurtosis and skewness of each run
% - gets an average overall skewness and kurtosis
% - finds the spread of these values (16th and 84th percentile)
% [1] Normal, Lognormal, Gamma, Weibull.
% We actually find that the random data predicts a mean and 'standard
% deviation' that can be compared directly to the data.

clear; clc; close all;
addpath("figures\kurtBias"); addpath("func/");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

%% Test 0-350 samples in steps of 5

% Specify parameters
N = [0.0456 0.0249];
L = [-3.3046 0.7849];
G = [2.4675 0.0185];
W = [0.0510 1.8469];

% Specify 'runs' r, i.e. how many times we repeat the randomly-distributed
% data calculation.
r = 5000;

% Sample Size s
s = 5:5:350;

for i=1:length(s)
    [n1(i),l1(i),sn1(i,:),sl1(i,:),g1(i),w1(i),pg1(i,:),pw1(i,:)] = testKurtosisBiasHelper(s(i),r,false,N,L,G,W);
    [mn1(i),ml1(i),ssn1(i,:),ssl1(i,:),g2(i),w2(i),pg2(i,:),pw2(i,:)] = testSkewnessBiasHelper(s(i),r,false,N,L,G,W);
end

%% Plot

ax1 = figure;
subplot(2,2,1)
plot(s,n1,Marker="+");
hold on
plot(s,sn1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(s,sn1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Normal");
subplot(2,2,2)
plot(s,l1);
hold on
plot(s,sl1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(s,sl1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Lognormal");
subplot(2,2,3)
plot(s,g1);
hold on
plot(s,pg1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(s,pg1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Gamma");
subplot(2,2,4)
plot(s,w1);
hold on
plot(s,pw1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(s,pw1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Weibull");
stitle = "Kurtosis Bias (" + sprintf('%d',r) + " r)";
sgtitle(stitle);

fName = "figures/kurtBias/__350_" + sprintf('%d',r) + "runs_C7.png";
exportgraphics(ax1,fName);

save("output\skku\kurtBiasC7.mat","n1","sn1","l1","sl1","g1","pg1","w1","pw1");

ax2 = figure;
subplot(2,2,1)
plot(s,mn1,Marker="+");
hold on
plot(s,ssn1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(s,ssn1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Normal");
subplot(2,2,2)
plot(s,ml1);
hold on
plot(s,ssl1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(s,ssl1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Lognormal");
subplot(2,2,3)
plot(s,g2);
hold on
plot(s,pg2(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(s,pg2(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Gamma");
subplot(2,2,4)
plot(s,w2);
hold on
plot(s,pw2(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(s,pw2(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Weibull");
stitle = "Skewness Bias (" + sprintf('%d',r) + " r)";
sgtitle(stitle);

fName = "figures/skewBias/__350_" + sprintf('%d',r) + "runs_C7.png";
exportgraphics(ax2,fName);

save("output\skku\skewBiasC7.mat","mn1","ssn1","ml1","ssl1","g2","pg2","w2","pw2");

clearvars -except r;