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
% deprecated, using stuff in isopycnalView.m now. But it might be worth
% running the 89-01 vs 01-21 comparison on the cruise-averaged night time
% code.

depthMeasurements = 129;
eulerianDepth = linspace(0,2*depthMeasurements,depthMeasurements);
lagrangianDepth = linspace(-128,128,depthMeasurements);
% 
% ksE = NaN(5,depthMeasurements);
% ksL = NaN(5,depthMeasurements);
% 
% % 12:329 = 1989 Oct -> 2021 Dec (complete - first year)
% % 12:196 = 1989 Oct -> 2007 Dec (BB's PhD)
% 
% % Eulerian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloro(i,12:329);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksE(:,i),~] = statsplot2(tmp,'noplot');
%     [h(i),p(i)] = lillietest(chloro(i,12:196));
%     [~,pln(i)] = lillietest(log(tmp));
%     [~,pwbl(i)] = lillietest(log(tmp),"Distr","ev");
% end
% clear tmp;
% 
% % Lagrangian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloroL(i,12:329);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksL(:,i),~] = statsplot2(tmp,'noplot');
% end
% 
% %% fi
% axa = figure;
% plot(ksE(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','KS','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(p,eulerianDepth,'DisplayName','Lillie');
% ylim([0 250]);
% set(gca,'YDir','reverse');
% legend();
% title('Goodness of Fit for a normal distribution (Eulerian)');
% exportgraphics(axa,'figures/ks-vs-lillie_norm_89-21.png');
% 
% axb = figure;
% plot(ksE(2,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','KS','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(pln,eulerianDepth,'DisplayName','Lillie');
% ylim([0 250]);
% set(gca,'YDir','reverse');
% legend();
% title('Goodness of Fit for a lognormal distribution (Eulerian)');
% exportgraphics(axb,'figures/ks-vs-lillie_logn_89-21.png');
% 
% axc = figure;
% plot(ksE(3,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','KS','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(pwbl,eulerianDepth,'DisplayName','Lillie');
% ylim([0 250]);
% set(gca,'YDir','reverse');
% legend();
% title('Goodness of Fit for a weibull distribution (Eulerian)');
% exportgraphics(axc,'figures/ks-vs-lillie_weib_89-21.png');
% 
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 9 15]);
% axd = figure;
% plot(p,eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(pln,eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(pwbl,eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% hold off
% ylim([0 250]);
% xlim([0 0.8]);
% set(gca,'YDir','reverse');
% legend();
% title('Lillie Test for Goodness of Fit');
% exportgraphics(axd,'figures/lillie_n_ln_wbl_89-21.png');
% 
% %% Plot the KS Statistic vs Depth for Five Distributions
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 15 15]);
% 
% ax1 = figure;
% % Eulerian
% subplot(1,2,1)
% plot(ksE(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksE(2,:),eulerianDepth,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksE(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksE(4,:),eulerianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% legend();
% ylim([0 250]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Eulerian KS-Test (89-21)');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksL(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksL(2,:),lagrangianDepth,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksL(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksL(4,:),lagrangianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% set(gca,'YDir','reverse');
% legend();
% ylim([-130 130]);
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Lagrangian KS-Test (89-21)');
% exportgraphics(ax1,'figures/ks-EulAndLag_1989_2021.png');
% 
% %% Is there a difference in KS if we exclude the first three years?
% 
% % Exclude September 1988 up to July 1991
% 
% ksE_3 = zeros(5,depthMeasurements);
% ksL_3 = zeros(5,depthMeasurements);
% newStartingPt = 28;
% 
% % Eulerian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloro(i,newStartingPt:end);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksE_3(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;
% 
% % Lagrangian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloroL(i,newStartingPt:end);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksL_3(:,i),~] = statsplot2(tmp,'noplot');
% end
% 
% %% Plot the KS Statistic vs Depth for Five Distributions (excl. 1988-1991)
% 
% ax2 = figure;
% % Eulerian
% subplot(1,2,1)
% plot(ksE_3(1,:),eulerianDepth,'Color',[0 0 0],'DisplayName','Normal');
% hold on
% plot(ksE_3(2,:),eulerianDepth,'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
% plot(ksE_3(3,:),eulerianDepth,'Color','red','DisplayName','Weibull');
% plot(ksE_3(4,:),eulerianDepth,'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksE_3(5,:),eulerianDepth,'Color','red','LineStyle',':','DisplayName','Exp');
% hold off
% legend();
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Eulerian KS-Test (excl. 1988 - 1991 Jul)');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksL_3(1,:),lagrangianDepth,'Color',[0 0 0],'DisplayName','Normal');
% hold on
% plot(ksL_3(2,:),lagrangianDepth,'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
% plot(ksL_3(3,:),lagrangianDepth,'Color','red','DisplayName','Weibull');
% plot(ksL_3(4,:),lagrangianDepth,'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksL_3(5,:),lagrangianDepth,'Color','red','LineStyle',':','DisplayName','Exp');
% hold off
% legend();
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Lagrangian KS-Test (excl. 1988 - 1991 Jul)');
% exportgraphics(ax2,'figures/ks-EulAndLag_1988_2021_excl3.png');
% %% ...and what about excluding the first five years?
% 
% % Exclude September 1988 up to September 1993
% 
% ksE_5 = zeros(5,depthMeasurements);
% ksL_5 = zeros(5,depthMeasurements);
% newStartingPt = 45;
% 
% % Eulerian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloro(i,newStartingPt:end);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksE_5(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;
% 
% % Lagrangian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloroL(i,newStartingPt:end);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksL_5(:,i),~] = statsplot2(tmp,'noplot');
% end
% 
% %% Plot the KS Statistic vs Depth for Five Distributions (excl. 1988-1991)
% 
% ax3 = figure;
% 
% % Eulerian
% subplot(1,2,1)
% plot(ksE_5(1,:),eulerianDepth,'Color',[0 0 0],'DisplayName','Normal');
% hold on
% plot(ksE_5(2,:),eulerianDepth,'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
% plot(ksE_5(3,:),eulerianDepth,'Color','red','DisplayName','Weibull');
% plot(ksE_5(4,:),eulerianDepth,'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksE_5(5,:),eulerianDepth,'Color','red','LineStyle',':','DisplayName','Exp');
% hold off
% legend();
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Eulerian KS-Test (excl. 1988 - 1991 Jul)');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksL_5(1,:),lagrangianDepth,'Color',[0 0 0],'DisplayName','Normal');
% hold on
% plot(ksL_5(2,:),lagrangianDepth,'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
% plot(ksL_5(3,:),lagrangianDepth,'Color','red','DisplayName','Weibull');
% plot(ksL_5(4,:),lagrangianDepth,'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksL_5(5,:),lagrangianDepth,'Color','red','LineStyle',':','DisplayName','Exp');
% hold off
% legend();
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Lagrangian KS-Test (excl. 1988 - 1993 Sep)');
% exportgraphics(ax3,'figures/ks-EulAndLag_1988_2021_excl5.png');
% 
% %% ...and just one year?
% 
% % Exclude September 1988 up to October 1989
% 
% ksE_1 = zeros(5,depthMeasurements);
% ksL_1 = zeros(5,depthMeasurements);
% newStartingPt = 12;
% 
% % Eulerian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloro(i,newStartingPt:end);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksE_1(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;
% 
% % Lagrangian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloroL(i,newStartingPt:end);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksL_1(:,i),~] = statsplot2(tmp,'noplot');
% end
% 
% %% Plot the KS Statistic vs Depth for Five Distributions (excl. 1988-1989)
% 
% ax4 = figure;
% 
% % Eulerian
% subplot(1,2,1)
% plot(ksE_1(1,:),eulerianDepth,'Color',[0 0 0],'DisplayName','Normal');
% hold on
% plot(ksE_1(2,:),eulerianDepth,'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
% plot(ksE_1(3,:),eulerianDepth,'Color','red','DisplayName','Weibull');
% plot(ksE_1(4,:),eulerianDepth,'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksE_1(5,:),eulerianDepth,'Color','red','LineStyle',':','DisplayName','Exp');
% hold off
% set(gca,'YDir','reverse');
% legend();
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Eulerian KS-Test (excl. 1988 - 1989 Oct)');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksL_1(1,:),lagrangianDepth,'Color',[0 0 0],'DisplayName','Normal');
% hold on
% plot(ksL_1(2,:),lagrangianDepth,'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
% plot(ksL_1(3,:),lagrangianDepth,'Color','red','DisplayName','Weibull');
% plot(ksL_1(4,:),lagrangianDepth,'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksL_1(5,:),lagrangianDepth,'Color','red','LineStyle',':','DisplayName','Exp');
% hold off
% set(gca,'YDir','reverse');
% legend();
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Lagrangian KS-Test (excl. 1988 - 1989 Oct)');
% exportgraphics(ax4,'figures/ks-EulAndLag_1988_2021_excl1.png');
% 
% %% Comparison of FOUR Possibilities
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 27 15]);
% 
% % Eulerian
% ax4 = figure;
% subplot(1,4,1)
% plot(ksE_1(1,:)-ksE(1,:),eulerianDepth,'DisplayName','89 Oct -');
% hold on
% plot(ksE_3(1,:)-ksE(1,:),eulerianDepth,'DisplayName','91 Jul -');
% plot(ksE_5(1,:)-ksE(1,:),eulerianDepth,'LineStyle','--','DisplayName','93 Feb -');
% hold off
% set(gca,'YDir','reverse');
% title('Normal');
% 
% subplot(1,4,2)
% plot(ksE_1(2,:)-ksE(2,:),eulerianDepth,'DisplayName','89 Oct -');
% hold on
% plot(ksE_3(2,:)-ksE(2,:),eulerianDepth,'DisplayName','91 Jul -');
% plot(ksE_5(2,:)-ksE(2,:),eulerianDepth,'LineStyle','--','DisplayName','93 Feb -');
% hold off
% set(gca,'YDir','reverse');
% title('Lognormal');
% 
% subplot(1,4,3)
% plot(ksE_1(3,:)-ksE(3,:),eulerianDepth,'DisplayName','89 Oct -');hold on
% plot(ksE_3(3,:)-ksE(3,:),eulerianDepth,'DisplayName','91 Jul -');
% plot(ksE_5(3,:)-ksE(3,:),eulerianDepth,'LineStyle','--','DisplayName','93 Feb -');hold off
% set(gca,'YDir','reverse');
% title('Weibull');
% 
% subplot(1,4,4)
% plot(ksE_1(4,:)-ksE(4,:),eulerianDepth,'DisplayName','89 Oct -');
% hold on
% plot(ksE_3(4,:)-ksE(4,:),eulerianDepth,'DisplayName','91 Jul -');
% plot(ksE_5(4,:)-ksE(4,:),eulerianDepth,'LineStyle','--','DisplayName','93 Feb -');
% hold off
% set(gca,'YDir','reverse');
% legend('Location','best');
% title('Gamma');
% 
% sgtitle('Change in KS Statistic based on excluding first one/three/five years');
% exportgraphics(ax4,'figures/ks-comparisonYearExclusions.png');
% 
% %% Compare KS for two eras: 1988 - 2001 and 2001 - 2021
% 
% % We compare these subsets because they represent different fluorometers.
% 
% ksE_era1 = zeros(5,depthMeasurements);
% ksL_era1 = zeros(5,depthMeasurements);
% ksE_era2 = zeros(5,depthMeasurements);
% ksL_era2 = zeros(5,depthMeasurements);
% 
% % 1988-2001 Eulerian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloro(i,1:128);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksE_era1(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;
% 
% % 1988-2001 Lagrangian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloroL(i,1:128);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksL_era1(:,i),~] = statsplot2(tmp,'noplot');
% end
% 
% % 2001-2021 Eulerian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloro(i,129:end);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksE_era2(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;
% 
% % 2001-2021 Lagrangian
% for i = 1:depthMeasurements
%     disp(i);
%     tmp = chloroL(i,129:end);
%     tmp(isnan(tmp) | tmp<=0) = [];
%     [~,ksLii(:,i),~] = statsplot2(tmp,'noplot');
% end
% 
% %%
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [7 2 16 18]);
% 
% % Eulerian
% ax5 = figure;
% subplot(1,2,1)
% plot(ksE_era1(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksE_era1(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksE_era1(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksE_era1(4,:),eulerianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% set(gca,'YDir','reverse');
% legend('Location','best');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Eulerian');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksL_era1(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksL_era1(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksL_era1(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksL_era1(4,:),lagrangianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Lagrangian');
% sgtitle('KS Test (1988-2001)');
% exportgraphics(ax5,'figures/ks_88_01.png');
% 
% ax6 = figure;
% % Eulerian
% subplot(1,2,1)
% plot(ksE_era2(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksE_era2(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksE_era2(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksE_era2(4,:),eulerianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% set(gca,'YDir','reverse');
% legend('Location','best');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Eulerian');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksLii(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksLii(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksLii(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksLii(4,:),lagrangianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('Depth [m]');
% title('Lagrangian');
% sgtitle('KS Test (2001-2021)')
% exportgraphics(ax6,'figures/ks_01_21.png');
% 
% %% Plot difference of 1988-2001 and 2001-2021 p-values
% 
% ax7 = figure;
% subplot(1,2,1)
% plot(ksE_era2(1,:)-ksE_era1(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksE_era2(2,:)-ksE_era1(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksE_era2(3,:)-ksE_era1(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksE_era2(4,:)-ksE_era1(4,:),eulerianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% set(gca,'YDir','reverse');
% legend('Location','best');
% title('Eulerian');
% subplot(1,2,2)
% plot(ksLii(1,:)-ksL_era1(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksLii(2,:)-ksL_era1(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksLii(3,:)-ksL_era1(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksLii(4,:)-ksL_era1(4,:),lagrangianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% set(gca,'YDir','reverse');
% title('Lagrangian');
% sgtitle('Change in p-value from 1988-2001 to 2001-2021');
% exportgraphics(ax7,'figures/ks_01_21_comparison.png');

%% Seasonal KS
% Need to apply this night time -> done in nightTimeKs.m

% Give numerical ID to each month
timeByMonth = discretize(month(time),12);

% Find ID for each season in the timeseries
winter = [];
spring = [];
summer = [];
autumn = [];

% 12:329 = 1989 Oct -> 2021 Dec (complete - first year)
% 12:196 = 1989 Oct -> 2007 Dec (BB's PhD)
for i=12:329
    if timeByMonth(i) <= 3
        winter = [winter i];
    elseif timeByMonth(i) >= 4 && timeByMonth(i) <= 6
        spring = [spring i];
    elseif timeByMonth(i) >= 7 && timeByMonth(i) <= 9
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

%% Seasonal KS: Eulerian Figures

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 25 18]);

ax8 = figure;
% WINTER
subplot(1,4,1)
plot(ksE_winter(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.2,'MarkerSize',4);
hold on
plot(ksE_winter(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.2,'MarkerSize',4);
plot(ksE_winter(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE_winter(4,:),eulerianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
legend('Location','best');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Winter (JFM)');

% SPRING
subplot(1,4,2)
plot(ksE_spring(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.2,'MarkerSize',4);
hold on
plot(ksE_spring(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.2,'MarkerSize',4);
plot(ksE_spring(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE_spring(4,:),eulerianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Spring (AMJ)');

% SUMMER
subplot(1,4,3)
plot(ksE_summer(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.2,'MarkerSize',4);
hold on
plot(ksE_summer(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.2,'MarkerSize',4);
plot(ksE_summer(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE_summer(4,:),eulerianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Summer (JAS)');

% AUTUMN
subplot(1,4,4)
plot(ksE_autumn(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksE_autumn(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksE_autumn(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE_autumn(4,:),eulerianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Autumn (OND)');

sgtitle('Kolmogorov-Smirnov Test: Seasonal, Eulerian (1989-2021)');
exportgraphics(ax8,'figures/ks_seasonal_eulerian_89-21.png');

%% Seasonal KS: Lagrangian Figures

ax9 = figure;
% WINTER
subplot(1,4,1)
plot(ksL_winter(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL_winter(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL_winter(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL_winter(4,:),lagrangianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
legend('Location','best');
set(gca, 'YDir','reverse');
ylim([-120 120]);
xlabel('p-value');
ylabel('Depth [m]');
title('Winter (JFM)');

% SPRING
subplot(1,4,2)
plot(ksL_spring(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL_spring(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL_spring(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL_spring(4,:),lagrangianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold 
set(gca, 'YDir','reverse');
ylim([-120 120]);
xlabel('p-value');
ylabel('Depth [m]');
title('Spring (AMJ)');

% SUMMER
subplot(1,4,3)
plot(ksL_summer(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL_summer(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL_summer(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL_summer(4,:),lagrangianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([-120 120]);
xlabel('p-value');
ylabel('Depth [m]');
title('Summer (JAS)');

% AUTUMN
subplot(1,4,4)
plot(ksL_autumn(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL_autumn(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL_autumn(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL_autumn(4,:),lagrangianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([-120 120]);
xlabel('p-value');
ylabel('Depth [m]');
title('Autumn (OND)');

sgtitle('Kolmogorov-Smirnov Test: Seasonal, Lagrangian (1989-2021)');
exportgraphics(ax9,'figures/ks_seasonal_lagrangian_89-21.png');

%% KS Test: stuff Benedetto just found!

% Eulerian
eulerianData = load('datafiles\hot_1773.mat');
fluoE = eulerianData.FLS_hot_1773(2:end,2:end);
pE = eulerianData.FLS_hot_1773(2:end,1);
ks_allcastE = nan(5,130);
for i = 1:130
    disp(i);
    tmp = fluoE(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    if length(tmp) > 2 
        [~,ks_allcastE(:,i),~] = statsplot2(tmp,'noplot');
    end
end

% Lagrangian
lagrangianData = load('datafiles\lagr_1758.mat');
fluoL = lagrangianData.FLS_lagr_1758(2:end,2:end);
pL = lagrangianData.FLS_lagr_1758(2:end,1);
ks_allcastL = nan(5,221);
for i = 1:221
    disp(i);
    tmp = fluoL(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    if length(tmp) > 2 
        [~,ks_allcastL(:,i),~] = statsplot2(tmp,'noplot');
    end
end

%% Plot Lagrangian of 1758 profiles

ax10 = figure;
subplot(1,2,1)
plot(ks_allcastE(1,:),pE,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4)
hold on
plot(ks_allcastE(2,:),pE,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks_allcastE(3,:),pE,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks_allcastE(4,:),pE,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 250]);
set(gca,'YDir','reverse');

title('Eulerian')

subplot(1,2,2)
plot(ks_allcastL(1,:),pL,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4)
hold on
plot(ks_allcastL(2,:),pL,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks_allcastL(3,:),pL,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks_allcastL(4,:),pL,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([-125 125]);
set(gca,'YDir','reverse');
title('Lagrangian');

ax10.Position = [3 3 20 15];
sgtitle('Yearly Kolmogorov-Smirnov Test');
exportgraphics(ax10,'figures/ks_89-07_allCasts.png');

%%
clear i axa axb axc axd ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10;