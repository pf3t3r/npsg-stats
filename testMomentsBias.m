clear; clc; close all;
addpath("figures\kurtBias");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

% This script will evaluate the performance of random data of varying
% sample size versus kurtosis and skewness. We aim to answer the following
% question: is there a bias for high kurtosis at low sample sizes? And is
% there a similar bias for skewness?

% Specify how many random arrays will be created in function.
runs = 100;

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
scatter(16,n1+sn1,'k',"+"); scatter(16,n1-sn1,'k',"+");
scatter(32,n2); scatter(32,n2+sn2,'k',"+"); scatter(32,n2-sn2,'k',"+");
scatter(64,n3); scatter(64,n3+sn3,'k',"+"); scatter(64,n3-sn3,'k',"+");
scatter(128,n4); scatter(128,n4+sn4,'k',"+"); scatter(128,n4-sn4,'k',"+");
scatter(256,n5); scatter(256,n5+sn5,'k',"+"); scatter(256,n5-sn5,'k',"+");
scatter(512,n6); scatter(512,n6+sn6,'k',"+"); scatter(512,n6-sn6,'k',"+");
scatter(1028,n7); scatter(1028,n7+sn7,'k',"+"); scatter(1028,n7-sn7,'k',"+");
scatter(2056,n8); scatter(2056,n8+sn8,'k',"+"); scatter(2056,n8-sn8,'k',"+");
scatter(4128,n9); scatter(4128,n9+sn9,'k',"+"); scatter(4128,n9-sn9,'k',"+");
scatter(8256,nX); scatter(8256,nX+snX,'k',"+"); scatter(8256,nX-snX,'k',"+");
scatter(16384,n11); scatter(16384,n11+sn11,'k',"+"); scatter(16384,n11-sn11,'k',"+");
scatter(32768,n12); scatter(32768,n12+sn12,'k',"+"); scatter(32768,n12-sn12,'k',"+");
scatter(65536,n13); scatter(65536,n13+sn13,'k',"+"); scatter(65536,n13-sn13,'k',"+");
scatter(131072,n14); scatter(131072,n14+sn14,'k',"+"); scatter(131072,n14-sn14,'k',"+");
hold off
grid on;
ylim([2 4]);
set(gca,'xscale','log');
title('Normal');

subplot(1,2,2)
scatter(16,l1);
hold on
scatter(16,l1+sl1,'k',"+"); scatter(16,l1-sl1,'k',"+");
scatter(32,l2); scatter(32,l2+sl2,'k',"+"); scatter(32,l2-sl2,'k',"+");
scatter(64,l3); scatter(64,l3+sl3,'k',"+"); scatter(64,l3-sl3,'k',"+");
scatter(128,l4); scatter(128,l4+sl4,'k',"+"); scatter(128,l4-sl4,'k',"+");
scatter(256,l5); scatter(256,l5+sl5,'k',"+"); scatter(256,l5-sl5,'k',"+");
scatter(512,l6); scatter(512,l6+sl6,'k',"+"); scatter(512,l6-sl6,'k',"+");
scatter(1028,l7); scatter(1028,l7+sl7,'k',"+"); scatter(1028,l7-sl7,'k',"+");
scatter(2056,l8); scatter(2056,l8+sl8,'k',"+"); scatter(2056,l8-sl8,'k',"+");
scatter(4128,l9); scatter(4128,l9+sl9,'k',"+"); scatter(4128,l9-sl9,'k',"+");
scatter(8256,lX); scatter(8256,lX+slX,'k',"+"); scatter(8256,lX-slX,'k',"+");
scatter(16384,l11); scatter(16384,l11+sl11,'k',"+"); scatter(16384,l11-sl11,'k',"+");
scatter(32768,l12); scatter(32768,l12+sl12,'k',"+"); scatter(32768,l12-sl12,'k',"+");
scatter(65536,l13); scatter(65536,l13+sl13,'k',"+"); scatter(65536,l13-sl13,'k',"+");
scatter(131072,l14); scatter(131072,l14+sl14,'k',"+"); scatter(131072,l14-sl14,'k',"+");
hold off
grid on;
ylim([1 7]);
set(gca,'xscale','log');
title('Lognormal');
stitle = "Kurtosis Bias (" + sprintf("%d",runs) + " runs)";
sgtitle(stitle);

fName = "figures/kurtBias/compOfBias_" + sprintf('%d',runs) + "runs.png";
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
scatter(16,mn1 + sn1,"k","+"); scatter(16,mn1 - sn1,"k","+");
scatter(32,mn2); scatter(32,mn2 + sn2,"k","+"); scatter(32,mn2 - sn2,"k","+");
scatter(64,mn3); scatter(64,mn3 + sn3,"k","+"); scatter(64,mn3 - sn3,"k","+");
scatter(128,mn4); scatter(128,mn4 + sn4,"k","+"); scatter(128,mn4 - sn4,"k","+");
scatter(256,mn5); scatter(256,mn5 + sn5,"k","+"); scatter(256,mn5 - sn5,"k","+");
scatter(512,mn6); scatter(512,mn6 + sn6,"k","+"); scatter(512,mn6 - sn6,"k","+");
scatter(1024,mn7); scatter(1024,mn7 + sn7,"k","+"); scatter(1024,mn7 - sn7,"k","+");
scatter(2048,mn8); scatter(2048,mn8 + sn8,"k","+"); scatter(2048,mn8 - sn8,"k","+");
scatter(4096,mn9); scatter(4096,mn9 + sn9,"k","+"); scatter(4096,mn9 - sn9,"k","+");
scatter(8192,mnX); scatter(8192,mnX + snX,"k","+"); scatter(8192,mnX - snX,"k","+");
scatter(16384,mn11); scatter(16384,mn11 + sn11,"k","+"); scatter(16384,mn11 - sn11,"k","+");
scatter(32768,mn12); scatter(32768,mn12 + sn12,"k","+"); scatter(32768,mn12 - sn12,"k","+");
scatter(65536,mn13); scatter(65536,mn13 + sn13,"k","+"); scatter(65536,mn13 - sn13,"k","+");
scatter(131072,mn14); scatter(131072,mn14 + sn14,"k","+"); scatter(131072,mn14 - sn14,"k","+");
hold off
grid on; ylim([-0.6 0.6]);
set(gca,'xscale','log');
title('Normal');

subplot(1,2,2)
scatter(16,ml1);
hold on
scatter(16,ml1 + sl1,"k","+"); scatter(16,ml1 - sl1,"k","+");
scatter(32,ml2); scatter(32,ml2 + sl2,"k","+"); scatter(32,ml2 - sl2,"k","+");
scatter(64,ml3); scatter(64,ml3 + sl3,"k","+"); scatter(64,ml3 - sl3,"k","+");
scatter(128,ml4); scatter(128,ml4 + sl4,"k","+"); scatter(128,ml4 - sl4,"k","+");
scatter(256,ml5); scatter(256,ml5 + sl5,"k","+"); scatter(256,ml5 - sl5,"k","+");
scatter(512,ml6); scatter(512,ml6 + sl6,"k","+"); scatter(512,ml6 - sl6,"k","+");
scatter(1024,ml7); scatter(1024,ml7 + sl7,"k","+"); scatter(1024,ml7 - sl7,"k","+");
scatter(2048,ml8); scatter(2048,ml8 + sl8,"k","+"); scatter(2048,ml8 - sl8,"k","+");
scatter(4096,ml9); scatter(4096,ml9 + sl9,"k","+"); scatter(4096,ml9 - sl9,"k","+");
scatter(8192,mlX); scatter(8192,mlX + slX,"k","+"); scatter(8192,mlX - slX,"k","+");
scatter(16384,ml11); scatter(16384,ml11 + sl11,"k","+"); scatter(16384,ml11 - sl11,"k","+");
scatter(32768,ml12); scatter(32768,ml12 + sl12,"k","+"); scatter(32768,ml12 - sl12,"k","+");
scatter(65536,ml13); scatter(65536,ml13 + sl13,"k","+"); scatter(65536,ml13 - sl13,"k","+");
scatter(131072,ml14); scatter(131072,ml14 + sl14,"k","+"); scatter(131072,ml14 - sl14,"k","+");
hold off
grid on; ylim([-0.2 1.6]);
set(gca,'xscale','log');
title('Lognormal');

stitle = "Skewness Bias (" + sprintf("%d",runs) + " runs)";
sgtitle(stitle);
fName = "figures/skewBias/compOfBias_" + sprintf('%d',runs) + "runs.png";
exportgraphics(ax,fName);

clearvars -except runs;

%% Test 0-250 samples in steps of 5

tmp = 5:5:350;
for i=1:length(tmp)
    [n1(i),l1(i),sn1(i),sl1(i)] = testKurtosisBiasHelper(tmp(i),runs);
    [mn1(i),ml1(i),ssn1(i),ssl1(i)] = testSkewnessBiasHelper(tmp(i),runs);
end

%%
ax1 = figure;
subplot(1,2,1)
plot(tmp,n1,Marker="+");
hold on
plot(tmp,n1+sn1,Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,n1-sn1,Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
subplot(1,2,2)
plot(tmp,l1);
hold on
plot(tmp,l1+sl1,Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,l1-sl1,Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
stitle = "Kurtosis Bias (" + sprintf('%d',runs) + " runs)";
sgtitle(stitle);

fName = "figures/kurtBias/_350_" + sprintf('%d',runs) + "runs.png";
exportgraphics(ax1,fName);

ax2 = figure;
subplot(1,2,1)
plot(tmp,mn1,Marker="+");
hold on
plot(tmp,mn1+ssn1,Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,mn1-ssn1,Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
subplot(1,2,2)
plot(tmp,ml1);
hold on
plot(tmp,ml1+ssl1,Marker=".",Color=[0.6 0.6 0.6]);
plot(tmp,ml1-ssl1,Marker=".",Color=[0.6 0.6 0.6]);
hold off
grid on;
stitle = "Skewness Bias (" + sprintf('%d',runs) + " runs)";
sgtitle(stitle);

fName = "figures/skewBias/_350_" + sprintf('%d',runs) + "runs.png";
exportgraphics(ax2,fName);

clearvars -except runs;

%% Test Z-values
tmp = 10:10:500;
for i=1:length(tmp)
    [kn(i),kl(i),knSD(i),klSD(i),knSE(i),klSE(i)] = testKurtosisBiasHelper(tmp(i),runs);
    %[mn1(i),ml1(i),ssn1(i),ssl1(i),seN(i),seL(i)] = testSkewnessBiasHelper(tmp(i),runs);
end

% z-values
zn = kn./knSE;
zl = kl./klSE;

figure;
plot(tmp,kn);
hold on
plot(tmp,kl);
hold off