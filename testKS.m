%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 15 15]);
set(0,'defaultAxesFontSize',10);

%% Load Data
chloro = load("datafiles\chloro.mat"','chloro256').chloro256;
chloroL = load("datafiles\chloro.mat"','chloro256_lang').chloro256_lang;
time = load("datafiles\chloro.mat","time256").time256;

%% Calculate the Kolmogorov-Smirnov Statistic for Chloropigment

depthMeasurements = 129;

mleE = zeros(5,2,depthMeasurements);
ksE = zeros(5,depthMeasurements);
nllE = zeros(5,depthMeasurements);

mleL = zeros(5,2,depthMeasurements);
ksL = zeros(5,depthMeasurements);
nllL = zeros(5,depthMeasurements);

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleE(:,:,i),ksE(:,i),nllE(:,i)] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleL(:,:,i),ksL(:,i),nllL(:,i)] = statsplot2(tmp,'noplot');
end

%% Plot the KS Statistic vs Depth for Five Distributions

ax1 = figure;
% Eulerian
subplot(1,2,1)
plot(ksE(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksE(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian KS-Test');

% Lagrangian
subplot(1,2,2)
plot(ksL(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksL(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian KS-Test');
exportgraphics(ax1,'figures/ks-EulAndLag_1988_2021.png');

%% Is there a difference in KS if we exclude the first three years?

% Exclude September 1988 up to July 1991

mleE_3 = zeros(5,2,depthMeasurements);
ksE_3 = zeros(5,depthMeasurements);
nllE_3 = zeros(5,depthMeasurements);

mleL_3 = zeros(5,2,depthMeasurements);
ksL_3 = zeros(5,depthMeasurements);
nllL_3 = zeros(5,depthMeasurements);

newStartingPt = 28;

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,newStartingPt:end);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleE_3(:,:,i),ksE_3(:,i),nllE_3(:,i)] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,newStartingPt:end);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleL_3(:,:,i),ksL_3(:,i),nllL_3(:,i)] = statsplot2(tmp,'noplot');
end

%% Plot the KS Statistic vs Depth for Five Distributions (excl. 1988-1991)

ax2 = figure;
% Eulerian
subplot(1,2,1)
plot(ksE_3(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE_3(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE_3(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE_3(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksE_3(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian KS-Test (excl. 1988 - 1991 Jul)');

% Lagrangian
subplot(1,2,2)
plot(ksL_3(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL_3(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL_3(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL_3(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksL_3(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian KS-Test (excl. 1988 - 1991 Jul)');
exportgraphics(ax2,'figures/ks-EulAndLag_1988_2021_excl3.png');
%% ...and what about excluding the first five years?

% Exclude September 1988 up to September 1993

mleE_5 = zeros(5,2,depthMeasurements);
ksE_5 = zeros(5,depthMeasurements);
nllE_5 = zeros(5,depthMeasurements);

mleL_5 = zeros(5,2,depthMeasurements);
ksL_5 = zeros(5,depthMeasurements);
nllL_5 = zeros(5,depthMeasurements);

newStartingPt = 45;

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,newStartingPt:end);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleE_5(:,:,i),ksE_5(:,i),nllE_5(:,i)] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,newStartingPt:end);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleL_5(:,:,i),ksL_5(:,i),nllL_5(:,i)] = statsplot2(tmp,'noplot');
end

%% Plot the KS Statistic vs Depth for Five Distributions (excl. 1988-1991)

ax3 = figure;

% Eulerian
subplot(1,2,1)
plot(ksE_5(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE_5(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE_5(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE_5(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksE_5(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian KS-Test (excl. 1988 - 1991 Jul)');

% Lagrangian
subplot(1,2,2)
plot(ksL_5(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL_5(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL_5(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL_5(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksL_5(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian KS-Test (excl. 1988 - 1993 Sep)');
exportgraphics(ax3,'figures/ks-EulAndLag_1988_2021_excl5.png');

%% ...and just one year?

% Exclude September 1988 up to October 1989

mleE_1 = zeros(5,2,depthMeasurements);
ksE_1 = zeros(5,depthMeasurements);
nllE_1 = zeros(5,depthMeasurements);

mleL_1 = zeros(5,2,depthMeasurements);
ksL_1 = zeros(5,depthMeasurements);
nllL_1 = zeros(5,depthMeasurements);

newStartingPt = 12;

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,newStartingPt:end);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleE_1(:,:,i),ksE_1(:,i),nllE_1(:,i)] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,newStartingPt:end);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleL_1(:,:,i),ksL_1(:,i),nllL_1(:,i)] = statsplot2(tmp,'noplot');
end

%% Plot the KS Statistic vs Depth for Five Distributions (excl. 1988-1989)

ax4 = figure;

% Eulerian
subplot(1,2,1)
plot(ksE_1(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE_1(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE_1(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE_1(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksE_1(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian KS-Test (excl. 1988 - 1989 Oct)');

% Lagrangian
subplot(1,2,2)
plot(ksL_1(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL_1(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL_1(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL_1(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksL_1(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian KS-Test (excl. 1988 - 1989 Oct)');
exportgraphics(ax4,'figures/ks-EulAndLag_1988_2021_excl1.png');

%% Comparison of FOUR Possibilities
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 27 15]);

% Eulerian
ax4 = figure;
subplot(1,5,1)
plot(ksE_1(1,:)-ksE(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','89 Oct -');
hold on
plot(ksE_3(1,:)-ksE(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','91 Jul -');
plot(ksE_5(1,:)-ksE(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'LineStyle','--','DisplayName','93 Feb -');
hold off
legend('Location','best');
title('Normal');

subplot(1,5,2)
plot(ksE_1(2,:)-ksE(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','89 Oct -');
hold on
plot(ksE_3(2,:)-ksE(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','91 Jul -');
plot(ksE_5(2,:)-ksE(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'LineStyle','--','DisplayName','93 Feb -');
hold off
legend('Location','best');
title('Lognormal');

subplot(1,5,3)
plot(ksE_1(3,:)-ksE(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','89 Oct -');hold on
plot(ksE_3(3,:)-ksE(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','91 Jul -');
plot(ksE_5(3,:)-ksE(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'LineStyle','--','DisplayName','93 Feb -');hold off
legend('Location','best');
title('Weibull');

subplot(1,5,4)
plot(ksE_1(4,:)-ksE(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','89 Oct -');
hold on
plot(ksE_3(4,:)-ksE(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','91 Jul -');
plot(ksE_5(4,:)-ksE(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'LineStyle','--','DisplayName','93 Feb -');
hold off
legend('Location','best');
title('Gamma');

subplot(1,5,5)
plot(ksE_1(5,:)-ksE(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','89 Oct -');
hold on
plot(ksE_3(5,:)-ksE(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'DisplayName','91 Jul -');
plot(ksE_5(5,:)-ksE(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'LineStyle','--','DisplayName','93 Feb -');
hold off
legend('Location','best');
title('Exp');

sgtitle('Change in KS Statistic based on excluding first 3 or 5 years of data');
exportgraphics(ax4,'figures/ks-comparisonYearExclusions.png');


%%
clear mleE mleL i ax1 ax2 ax3 ax4;