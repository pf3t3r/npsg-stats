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
chloroL = load("datafiles\chloro.mat"','chloro256l').chloro256l;
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

sgtitle('Change in KS Statistic based on excluding first one/three/five years');
exportgraphics(ax4,'figures/ks-comparisonYearExclusions.png');

%% Compare KS for two subsets: 1988 - 2001 and 2001 - 2021

% We compare these subsets because they represent different fluorometers.

mleEi = zeros(5,2,depthMeasurements);
ksEi = zeros(5,depthMeasurements);
nllEi = zeros(5,depthMeasurements);
mleLi = zeros(5,2,depthMeasurements);
ksLi = zeros(5,depthMeasurements);
nllLi = zeros(5,depthMeasurements);
mleEii = zeros(5,2,depthMeasurements);
ksEii = zeros(5,depthMeasurements);
nllEii = zeros(5,depthMeasurements);
mleLii = zeros(5,2,depthMeasurements);
ksLii = zeros(5,depthMeasurements);
nllLii = zeros(5,depthMeasurements);

% 1988-2001 Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,1:128);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleEi(:,:,i),ksEi(:,i),nllEi(:,i)] = statsplot2(tmp,'noplot');
end
clear tmp;

% 1988-2001 Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,1:128);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleLi(:,:,i),ksLi(:,i),nllLi(:,i)] = statsplot2(tmp,'noplot');
end

% 2001-2021 Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,129:end);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleEii(:,:,i),ksEii(:,i),nllEii(:,i)] = statsplot2(tmp,'noplot');
end
clear tmp;

% 2001-2021 Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,129:end);
    tmp(isnan(tmp) | tmp<=0) = [];
    [mleLii(:,:,i),ksLii(:,i),nllLii(:,i)] = statsplot2(tmp,'noplot');
end

% Eulerian
ax5 = figure;
subplot(1,2,1)
plot(ksEi(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksEi(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksEi(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksEi(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksEi(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian');

% Lagrangian
subplot(1,2,2)
plot(ksLi(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksLi(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksLi(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksLi(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksLi(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian');
sgtitle('KS Test (1988-2001)');
exportgraphics(ax5,'figures/ks_88_01.png');

ax6 = figure;
% Eulerian
subplot(1,2,1)
plot(ksEii(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksEii(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksEii(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksEii(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksEii(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian');

% Lagrangian
subplot(1,2,2)
plot(ksLii(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksLii(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksLii(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksLii(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksLii(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian');
sgtitle('KS Test (2001-2021)')
exportgraphics(ax6,'figures/ks_01_21.png');

%% Plot difference of 1988-2001 and 2001-2021 p-values

ax7 = figure;
subplot(1,2,1)
plot(ksEii(1,:)-ksEi(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksEii(2,:)-ksEi(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksEii(3,:)-ksEi(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksEii(4,:)-ksEi(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksEii(5,:)-ksEi(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','blue','DisplayName','Exp');
hold off
legend();
title('Eulerian');
subplot(1,2,2)
plot(ksLii(1,:)-ksLi(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksLii(2,:)-ksLi(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksLii(3,:)-ksLi(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksLii(4,:)-ksLi(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(ksLii(5,:)-ksLi(5,:),linspace(128,-128,depthMeasurements),'Color','blue','DisplayName','Exp');
hold off
legend();
title('Lagrangian');
sgtitle('Change in p-value from 1988-2001 to 2001-2021');
exportgraphics(ax7,'figures/ks_01_21_comparison.png');

%% Seasonal KS

% chloroSpring = 
% yourtimetable.Season = discretize(mod(month(test1.Time), 12), 0:3:12, 'categorical', {'winter', 'spring', 'summer', 'autumn'});
%the mod is to bring december as first month

% Give numerical ID to each month
timeByMonth = discretize(month(time),12);

% Find ID for each season in the timeseries
winter = [];
spring = [];
summer = [];
autumn = [];

for i=1:329
    if timeByMonth(i) == 12 || timeByMonth(i) <= 2
        winter = [winter i];
    elseif timeByMonth(i) >= 3 && timeByMonth(i) <= 5
        spring = [spring i];
    elseif timeByMonth(i) >= 6 && timeByMonth(i) <= 8
        summer = [summer i];
    else
        autumn = [autumn i];
    end
end

%% Apply KS to seasonal data

% Try Eulerian first.
ksE_winter = zeros(5,depthMeasurements);
ksE_spring = zeros(5,depthMeasurements);
ksE_summer = zeros(5,depthMeasurements);
ksE_autumn = zeros(5,depthMeasurements);

% Try Lagrangian later.
ksL_winter = zeros(5,depthMeasurements);
ksL_spring = zeros(5,depthMeasurements);
ksL_summer = zeros(5,depthMeasurements);
ksL_autumn = zeros(5,depthMeasurements);

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,winter);
    tmp(isnan(tmp) | tmp<=0) = [];
    [~,ksE_winter(:,i),~] = statsplot2(tmp,'noplot');
    tmp2 = chloro(i,spring);
    tmp2(isnan(tmp2) | tmp2<=0) = [];
    [~,ksE_spring(:,i),~] = statsplot2(tmp2,'noplot');
    tmp3 = chloro(i,summer);
    tmp3(isnan(tmp3) | tmp3<=0) = [];
    [~,ksE_summer(:,i),~] = statsplot2(tmp3,'noplot');
    tmp4 = chloro(i,autumn);
    tmp4(isnan(tmp4) | tmp4<=0) = [];
    [~,ksE_autumn(:,i),~] = statsplot2(tmp4,'noplot');
end
clear tmp;


% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,winter);
    tmp(isnan(tmp) | tmp<=0) = [];
    [~,ksL_winter(:,i),~] = statsplot2(tmp,'noplot');
    tmp2 = chloroL(i,spring);
    tmp2(isnan(tmp2) | tmp2<=0) = [];
    [~,ksL_spring(:,i),~] = statsplot2(tmp2,'noplot');
    tmp3 = chloroL(i,summer);
    tmp3(isnan(tmp3) | tmp3<=0) = [];
    [~,ksL_summer(:,i),~] = statsplot2(tmp3,'noplot');
    tmp4 = chloroL(i,autumn);
    tmp4(isnan(tmp4) | tmp4<=0) = [];
    [~,ksL_autumn(:,i),~] = statsplot2(tmp4,'noplot');
end

%% Plot the seasonal KS for Eulerian
ax8 = figure;
% Eulerian

% WINTER
subplot(1,4,1)
plot(ksE_winter(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE_winter(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE_winter(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE_winter(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksE_winter(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Winter (DJF)');

% SPRING
subplot(1,4,2)
plot(ksE_spring(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE_spring(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE_spring(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE_spring(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksE_spring(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Spring (MAM)');

% SUMMER
subplot(1,4,3)
plot(ksE_summer(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE_summer(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE_summer(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE_summer(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksE_summer(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Summer (JJA)');

% AUTUMN
subplot(1,4,4)
plot(ksE_autumn(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksE_autumn(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksE_autumn(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksE_autumn(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksE_autumn(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Autumn (SON)');

sgtitle('Kolmogorov-Smirnov Test: Seasonal, Eulerian (1988-2021)');
exportgraphics(ax8,'figures/ks_seasonal_eulerian.png');

%% Plot the seasonal KS for Lagrangian
ax9 = figure;
% Lagrangian

% WINTER
subplot(1,4,1)
plot(ksL_winter(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL_winter(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL_winter(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL_winter(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksL_winter(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Winter (DJF)');

% SPRING
subplot(1,4,2)
plot(ksL_spring(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL_spring(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL_spring(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL_spring(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksE_spring(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Spring (MAM)');

% SUMMER
subplot(1,4,3)
plot(ksL_summer(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL_summer(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL_summer(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL_summer(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksL_summer(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Summer (JJA)');

% AUTUMN
subplot(1,4,4)
plot(ksL_autumn(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksL_autumn(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksL_autumn(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(ksL_autumn(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksL_autumn(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Autumn (SON)');

sgtitle('Kolmogorov-Smirnov Test: Seasonal, Lagrangian (1988-2021)');
exportgraphics(ax9,'figures/ks_seasonal_lagrangian.png');

%%
clear mleE mleL i ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9;