% Check relationship between p-value and sample size.
% This code checks whether there is any relationship between the p-value of
% three hypothesis tests [1] run on random lognormal data and the (i) mean* 
% and (ii) sample size [2] chosen for that data. It appears that there is 
% no significant dependence on mean or sample size. It may be worth
% rerunning the test with a more representative (of chl-a) parameter
% selection or with varying STD*.
% [1] K-S, Lil, A-D.
% [2] See variables a1 - a6.
% *For this test, the STD parameter was kept constant.

clc; close all; clear;

%% Initialise MEAN and STD

mu = [0.5 1 1.5 2 2.5 3];
sd = 0.4;

%% Initialise random arrays

% Adapt number initialised
a1 = 500; a2 = 1000; a3 = 1500; a4 = 2000; a5 = 2500; a6 = 3000;

tmp20 = []; tmp40 = []; tmp60 = []; tmp80 = []; tmp100 = []; tmp1000 = [];

tmp20 = nan(a1,50,6);
tmp40 = nan(a2,50,6);
tmp60 = nan(a3,50,6);
tmp80 = nan(a4,50,6);
tmp100 = nan(a5,50,6);
tmp1000 = nan(a6,50,6);

for j = 1:50
for i = 1:length(mu)
    t20 = lognrnd(mu(i),sd,a1,1);
    tmp20(:,j,i) = t20;
    t40 = lognrnd(mu(i),sd,a2,1);
    tmp40(:,j,i) = t40;
    t60 = lognrnd(mu(i),sd,a3,1);
    tmp60(:,j,i) = t60;
    t80 = lognrnd(mu(i),sd,a4,1);
    tmp80(:,j,i) = t80;
    t100 = lognrnd(mu(i),sd,a5,1);
    tmp100(:,j,i) = t100;
    t1000 = lognrnd(mu(i),sd,a6,1);
    tmp1000(:,j,i) = t1000;
end
end 

%% K-S, Lil, A-D: Lognormal test of lognormal data

for i = 1:50
    [ksL20(i,:),lilL20(i,:),adL20(i,:)] = quickMleAndKs(tmp20(:,i,:),6);
    [ksL40(i,:),lilL40(i,:),adL40(i,:)] = quickMleAndKs(tmp40(:,i,:),6);
    [ksL60(i,:),lilL60(i,:),adL60(i,:)] = quickMleAndKs(tmp60(:,i,:),6);
    [ksL80(i,:),lilL80(i,:),adL80(i,:)] = quickMleAndKs(tmp80(:,i,:),6);
    [ksL100(i,:),lilL100(i,:),adL100(i,:)] = quickMleAndKs(tmp100(:,i,:),6);
    [ksL1000(i,:),lilL1000(i,:),adL1000(i,:)] = quickMleAndKs(tmp1000(:,i,:),6);
end

% Save (Optional)
% save("output\Prec\mu6.mat","ksL20","ksL40","ksL60","ksL80","ksL100","ksL1000");

%% K-S Lognormal Test of Lognormal Data

ax = figure;
plot(mu,mean(ksL20),'Color','#a6cee3','DisplayName',sprintf('n = %d', a1));
hold on
plot(mu,mean(ksL20)+std(ksL20),'Color','#a6cee3','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL20)-std(ksL20),'Color','#a6cee3','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL40),'Color','#1f78b4','DisplayName',sprintf('n = %d', a2));
plot(mu,mean(ksL40)+std(ksL40),'Color','#1f78b4','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL40)-std(ksL40),'Color','#1f78b4','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL60),'Color','#b2df8a','DisplayName',sprintf('n = %d', a3));
plot(mu,mean(ksL60)+std(ksL60),'Color','#b2df8a','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL60)-std(ksL60),'Color','#b2df8a','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL80),'Color','#33a02c','DisplayName',sprintf('n = %d', a4));
plot(mu,mean(ksL80)+std(ksL80),'Color','#33a02c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL80)-std(ksL80),'Color','#33a02c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL100),'Color','#fb9a99','DisplayName',sprintf('n = %d', a5));
plot(mu,mean(ksL100)+std(ksL100),'Color','#fb9a99','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL100)-std(ksL100),'Color','#fb9a99','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL1000),'Color','#e31a1c','DisplayName',sprintf('n = %d', a6));
plot(mu,mean(ksL1000)+std(ksL1000),'Color','#e31a1c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(ksL1000)-std(ksL1000),'Color','#e31a1c','LineStyle',':',HandleVisibility="off");
hold off
ylabel('p-value');
xlabel('mean');
legend(Location="best");
title('K-S Lognormal Test of Random Lognormal Data','Sensitivity Analysis (\sigma = 0.4)');
exportgraphics(ax,'figures/prec/runKs_xHigh.png'); clear ax;

%% Lilliefors Lognormal Test of Lognormal Data
% It appears as though the p-value is not returned for values over 0.5. Not
% sure why yet.

ax2 = figure;
plot(mu,mean(lilL20),'Color','#a6cee3','DisplayName',sprintf('n = %d', a1));
hold on
plot(mu,mean(lilL20)+std(lilL20),'Color','#a6cee3','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL20)-std(lilL20),'Color','#a6cee3','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL40),'Color','#1f78b4','DisplayName',sprintf('n = %d', a2));
plot(mu,mean(lilL40)+std(lilL40),'Color','#1f78b4','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL40)-std(lilL40),'Color','#1f78b4','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL60),'Color','#b2df8a','DisplayName',sprintf('n = %d', a3));
plot(mu,mean(lilL60)+std(lilL60),'Color','#b2df8a','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL60)-std(lilL60),'Color','#b2df8a','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL80),'Color','#33a02c','DisplayName',sprintf('n = %d', a4));
plot(mu,mean(lilL80)+std(lilL80),'Color','#33a02c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL80)-std(lilL80),'Color','#33a02c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL100),'Color','#fb9a99','DisplayName',sprintf('n = %d', a5));
plot(mu,mean(lilL100)+std(lilL100),'Color','#fb9a99','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL100)-std(lilL100),'Color','#fb9a99','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL1000),'Color','#e31a1c','DisplayName',sprintf('n = %d', a6));
plot(mu,mean(lilL1000)+std(lilL1000),'Color','#e31a1c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(lilL1000)-std(lilL1000),'Color','#e31a1c','LineStyle',':',HandleVisibility="off");
ylabel('p-value');
xlabel('mean');
legend(Location="best");
title('Lilliefors Lognormal Test of Random Lognormal Data','Sensitivity Analysis (\sigma = 0.4)');
exportgraphics(ax2,'figures/prec/runLil_xHigh.png'); clear ax2;

%% A-D Lognormal Test of Lognormal Data

ax3 = figure;
plot(mu,mean(adL20),'Color','#a6cee3','DisplayName',sprintf('n = %d', a1));
hold on
plot(mu,mean(adL20)+std(adL20),'Color','#a6cee3','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL20)-std(adL20),'Color','#a6cee3','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL40),'Color','#1f78b4','DisplayName',sprintf('n = %d', a2));
plot(mu,mean(adL40)+std(adL40),'Color','#1f78b4','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL40)-std(adL40),'Color','#1f78b4','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL60),'Color','#b2df8a','DisplayName',sprintf('n = %d', a3));
plot(mu,mean(adL60)+std(adL60),'Color','#b2df8a','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL60)-std(adL60),'Color','#b2df8a','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL80),'Color','#33a02c','DisplayName',sprintf('n = %d', a4));
plot(mu,mean(adL80)+std(adL80),'Color','#33a02c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL80)-std(adL80),'Color','#33a02c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL100),'Color','#fb9a99','DisplayName',sprintf('n = %d', a5));
plot(mu,mean(adL100)+std(adL100),'Color','#fb9a99','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL100)-std(adL100),'Color','#fb9a99','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL1000),'Color','#e31a1c','DisplayName',sprintf('n = %d', a6));
plot(mu,mean(adL1000)+std(adL1000),'Color','#e31a1c','LineStyle',':',HandleVisibility="off");
plot(mu,mean(adL1000)-std(adL1000),'Color','#e31a1c','LineStyle',':',HandleVisibility="off");
ylabel('p-value');
xlabel('mean');
legend(Location="best");
title('A-D Lognormal Test of Random Lognormal Data','Sensitivity Analysis (\sigma = 0.4)');
exportgraphics(ax3,'figures/prec/runAd_xHigh.png'); clear ax3;