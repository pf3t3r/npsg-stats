clear; clc; close all;
addpath("figures\kurtBias");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

%% TODO: Add Weibull and Gamma Analysis.

%% What does this script do?

% This script will evaluate the performance of random data of varying
% sample size versus kurtosis and skewness. We aim to answer the following
% question: is there a bias for high kurtosis at low sample sizes? And is
% there a similar bias for skewness?

% Specify how many random arrays will be created in function.
runs = 5000;

%% Kurtosis Bias

[n1,l1,sn1,sl1] = testKurtosisBiasHelper(16,runs);
[n2,l2,sn2,sl2] = testKurtosisBiasHelper(32,runs);
[n3,l3,sn3,sl3] = testKurtosisBiasHelper(64,runs);
[n4,l4,sn4,sl4] = testKurtosisBiasHelper(128,runs);
[n5,l5,sn5,sl5] = testKurtosisBiasHelper(256,runs);
[n6,l6,sn6,sl6] = testKurtosisBiasHelper(512,runs);
[n7,l7,sn7,sl7] = testKurtosisBiasHelper(1024,runs);
[n8,l8,sn8,sl8] = testKurtosisBiasHelper(2048,runs);
[n9,l9,sn9,sl9] = testKurtosisBiasHelper(4096,runs);
[nX,lX,snX,slX] = testKurtosisBiasHelper(8192,runs);
[n11,l11,sn11,sl11] = testKurtosisBiasHelper(16384,runs);
[n12,l12,sn12,sl12] = testKurtosisBiasHelper(32768,runs);
[n13,l13,sn13,sl13] = testKurtosisBiasHelper(65536,runs);
[n14,l14,sn14,sl14] = testKurtosisBiasHelper(131072,runs);

ax = figure;
subplot(1,2,1)
scatter(16,n1);
hold on
scatter(16,sn1(1),'k',"+"); scatter(16,sn1(2),'k',"+");
scatter(32,n2); scatter(32,sn2(1),'k',"+"); scatter(32,sn2(2),'k',"+");
scatter(64,n3); scatter(64,sn3(1),'k',"+"); scatter(64,sn3(2),'k',"+");
scatter(128,n4); scatter(128,sn4(1),'k',"+"); scatter(128,sn4(2),'k',"+");
scatter(256,n5); scatter(256,sn5(1),'k',"+"); scatter(256,sn5(2),'k',"+");
scatter(512,n6); scatter(512,sn6(1),'k',"+"); scatter(512,sn6(2),'k',"+");
scatter(1028,n7); scatter(1028,sn7(1),'k',"+"); scatter(1028,sn7(2),'k',"+");
scatter(2056,n8); scatter(2056,sn8(1),'k',"+"); scatter(2056,sn8(2),'k',"+");
scatter(4128,n9); scatter(4128,sn9(1),'k',"+"); scatter(4128,sn9(2),'k',"+");
scatter(8256,nX); scatter(8256,snX(1),'k',"+"); scatter(8256,snX(2),'k',"+");
scatter(16384,n11); scatter(16384,sn11(1),'k',"+"); scatter(16384,sn11(2),'k',"+");
scatter(32768,n12); scatter(32768,sn12(1),'k',"+"); scatter(32768,sn12(2),'k',"+");
scatter(65536,n13); scatter(65536,sn13(1),'k',"+"); scatter(65536,sn13(2),'k',"+");
scatter(131072,n14); scatter(131072,sn14(1),'k',"+"); scatter(131072,sn14(2),'k',"+");
hold off
grid on;
ylim([2 4]);
set(gca,'xscale','log');
title('Normal');

subplot(1,2,2)
scatter(16,l1);
hold on
scatter(16,sl1(1),'k',"+"); scatter(16,sl1(2),'k',"+");
scatter(32,l2); scatter(32,sl2(1),'k',"+"); scatter(32,sl2(2),'k',"+");
scatter(64,l3); scatter(64,sl3(1),'k',"+"); scatter(64,sl3(2),'k',"+");
scatter(128,l4); scatter(128,sl4(1),'k',"+"); scatter(128,sl4(2),'k',"+");
scatter(256,l5); scatter(256,sl5(1),'k',"+"); scatter(256,sl5(2),'k',"+");
scatter(512,l6); scatter(512,sl6(1),'k',"+"); scatter(512,sl6(2),'k',"+");
scatter(1028,l7); scatter(1028,sl7(1),'k',"+"); scatter(1028,sl7(2),'k',"+");
scatter(2056,l8); scatter(2056,sl8(1),'k',"+"); scatter(2056,sl8(2),'k',"+");
scatter(4128,l9); scatter(4128,sl9(1),'k',"+"); scatter(4128,sl9(2),'k',"+");
scatter(8256,lX); scatter(8256,slX(1),'k',"+"); scatter(8256,slX(2),'k',"+");
scatter(16384,l11); scatter(16384,sl11(1),'k',"+"); scatter(16384,sl11(2),'k',"+");
scatter(32768,l12); scatter(32768,sl12(1),'k',"+"); scatter(32768,sl12(2),'k',"+");
scatter(65536,l13); scatter(65536,sl13(1),'k',"+"); scatter(65536,sl13(2),'k',"+");
scatter(131072,l14); scatter(131072,sl14(1),'k',"+"); scatter(131072,sl14(2),'k',"+");
hold off
grid on;
ylim([1 7]);
set(gca,'xscale','log');
title('Lognormal');
stitle = "Kurtosis Bias (" + sprintf("%d",runs) + " runs)";
sgtitle(stitle);

fName = "figures/kurtBias/__" + sprintf('%d',runs) + "runs.png";
exportgraphics(ax,fName);

clearvars -except runs;

%% Skewness Bias

[mn1,ml1,sn1,sl1] = testSkewnessBiasHelper(16,runs);
[mn2,ml2,sn2,sl2] = testSkewnessBiasHelper(32,runs);
[mn3,ml3,sn3,sl3] = testSkewnessBiasHelper(64,runs);
[mn4,ml4,sn4,sl4] = testSkewnessBiasHelper(128,runs);
[mn5,ml5,sn5,sl5] = testSkewnessBiasHelper(256,runs);
[mn6,ml6,sn6,sl6] = testSkewnessBiasHelper(512,runs);
[mn7,ml7,sn7,sl7] = testSkewnessBiasHelper(1024,runs);
[mn8,ml8,sn8,sl8] = testSkewnessBiasHelper(2048,runs);
[mn9,ml9,sn9,sl9] = testSkewnessBiasHelper(4096,runs);
[mnX,mlX,snX,slX] = testSkewnessBiasHelper(8192,runs);
[mn11,ml11,sn11,sl11] = testSkewnessBiasHelper(16384,runs);
[mn12,ml12,sn12,sl12] = testSkewnessBiasHelper(32768,runs);
[mn13,ml13,sn13,sl13] = testSkewnessBiasHelper(65536,runs);
[mn14,ml14,sn14,sl14] = testSkewnessBiasHelper(131072,runs);

ax = figure;
subplot(1,2,1)
scatter(16,mn1);
hold on
scatter(16,sn1(1),"k","+"); scatter(16,sn1(2),"k","+");
scatter(32,mn2); scatter(32,sn2(1),"k","+"); scatter(32,sn2(2),"k","+");
scatter(64,mn3); scatter(64,sn3(1),"k","+"); scatter(64,sn3(2),"k","+");
scatter(128,mn4); scatter(128,sn4(1),"k","+"); scatter(128,sn4(2),"k","+");
scatter(256,mn5); scatter(256,sn5(1),"k","+"); scatter(256,sn5(2),"k","+");
scatter(512,mn6); scatter(512,sn6(1),"k","+"); scatter(512,sn6(2),"k","+");
scatter(1024,mn7); scatter(1024,sn7(1),"k","+"); scatter(1024,sn7(2),"k","+");
scatter(2048,mn8); scatter(2048,sn8(1),"k","+"); scatter(2048,sn8(2),"k","+");
scatter(4096,mn9); scatter(4096,sn9(1),"k","+"); scatter(4096,sn9(2),"k","+");
scatter(8192,mnX); scatter(8192,snX(1),"k","+"); scatter(8192,snX(2),"k","+");
scatter(16384,mn11); scatter(16384,sn11(1),"k","+"); scatter(16384,sn11(2),"k","+");
scatter(32768,mn12); scatter(32768,sn12(1),"k","+"); scatter(32768,sn12(2),"k","+");
scatter(65536,mn13); scatter(65536,sn13(1),"k","+"); scatter(65536,sn13(2),"k","+");
scatter(131072,mn14); scatter(131072,sn14(1),"k","+"); scatter(131072,sn14(2),"k","+");
hold off
grid on; ylim([-0.6 0.6]);
set(gca,'xscale','log');
title('Normal');

subplot(1,2,2)
scatter(16,ml1);
hold on
scatter(16,sl1(1),"k","+"); scatter(16,sl1(2),"k","+");
scatter(32,ml2); scatter(32,sl2(1),"k","+"); scatter(32,sl2(2),"k","+");
scatter(64,ml3); scatter(64,sl3(1),"k","+"); scatter(64,sl3(2),"k","+");
scatter(128,ml4); scatter(128,sl4(1),"k","+"); scatter(128,sl4(2),"k","+");
scatter(256,ml5); scatter(256,sl5(1),"k","+"); scatter(256,sl5(2),"k","+");
scatter(512,ml6); scatter(512,sl6(1),"k","+"); scatter(512,sl6(2),"k","+");
scatter(1024,ml7); scatter(1024,sl7(1),"k","+"); scatter(1024,sl7(2),"k","+");
scatter(2048,ml8); scatter(2048,sl8(1),"k","+"); scatter(2048,sl8(2),"k","+");
scatter(4096,ml9); scatter(4096,sl9(1),"k","+"); scatter(4096,sl9(2),"k","+");
scatter(8192,mlX); scatter(8192,slX(1),"k","+"); scatter(8192,slX(2),"k","+");
scatter(16384,ml11); scatter(16384,sl11(1),"k","+"); scatter(16384,sl11(2),"k","+");
scatter(32768,ml12); scatter(32768,sl12(1),"k","+"); scatter(32768,sl12(2),"k","+");
scatter(65536,ml13); scatter(65536,sl13(1),"k","+"); scatter(65536,sl13(2),"k","+");
scatter(131072,ml14); scatter(131072,sl14(1),"k","+"); scatter(131072,sl14(2),"k","+");
hold off
grid on; ylim([-0.2 1.6]);
set(gca,'xscale','log');
title('Lognormal');

stitle = "Skewness Bias (" + sprintf("%d",runs) + " runs)";
sgtitle(stitle);
fName = "figures/skewBias/__" + sprintf('%d',runs) + "runs.png";
exportgraphics(ax,fName);

clearvars -except runs;

%% Test 0-250 samples in steps of 5

tmp = 5:5:350;
for i=1:length(tmp)
    [n1(i),l1(i),sn1(i,:),sl1(i,:),g1(i),w1(i),pg1(i,:),pw1(i,:)] = testKurtosisBiasHelper(tmp(i),runs);
    [mn1(i),ml1(i),ssn1(i,:),ssl1(i,:),g2(i),w2(i),pg2(i,:),pw2(i,:)] = testSkewnessBiasHelper(tmp(i),runs);
end

%%
ax1 = figure;
subplot(2,2,1)
plot(tmp,n1,Marker="+");
hold on
plot(tmp,sn1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,sn1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Normal");
subplot(2,2,2)
plot(tmp,l1);
hold on
plot(tmp,sl1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,sl1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Lognormal");
subplot(2,2,3)
plot(tmp,g1);
hold on
plot(tmp,pg1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,pg1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Gamma");
subplot(2,2,4)
plot(tmp,w1);
hold on
plot(tmp,pw1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,pw1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Weibull");
stitle = "Kurtosis Bias (" + sprintf('%d',runs) + " runs)";
sgtitle(stitle);

fName = "figures/kurtBias/__350_" + sprintf('%d',runs) + "runs.png";
exportgraphics(ax1,fName);

save("output\skku\kurtBias.mat","n1","sn1","l1","sl1","g1","pg1","w1","pw1");

ax2 = figure;
subplot(2,2,1)
plot(tmp,mn1,Marker="+");
hold on
plot(tmp,ssn1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,ssn1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Normal");
subplot(2,2,2)
plot(tmp,ml1);
hold on
plot(tmp,ssl1(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,ssl1(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Lognormal");
subplot(2,2,3)
plot(tmp,g2);
hold on
plot(tmp,pg2(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,pg2(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Gamma");
subplot(2,2,4)
plot(tmp,w2);
hold on
plot(tmp,pw2(:,1),Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,pw2(:,2),Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
title("Weibull");
stitle = "Skewness Bias (" + sprintf('%d',runs) + " runs)";
sgtitle(stitle);

fName = "figures/skewBias/__350_" + sprintf('%d',runs) + "runs.png";
exportgraphics(ax2,fName);

save("output\skku\skewBias.mat","mn1","ssn1","ml1","ssl1","g2","pg2","w2","pw2");


clearvars -except runs;