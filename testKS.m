clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
% set(groot,'defaultAxesXGrid','on');
% set(groot,'defaultAxesYGrid','on');
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 28 12]);
set(0,'defaultAxesFontSize',10);

%% Load Test Data
chloro = load("datafiles\chloro.mat"','chloro256').chloro256;
chloroL = load("datafiles\chloro.mat"','chloro256_lang').chloro256_lang;
time = load("datafiles\chloro.mat","time256").time256;

%% Remove NaN or Zeros

tickyboi = [];
antiTick = [];
tickyboiL = [];
antiTickL = [];

for i = 1:329
    for j = 1:129
        if (~any(chloro(j,i)))
            tickyboi = [tickyboi i];
            break;
        else
            antiTick = [antiTick i];
            break;
        end
    end
end


for i = 1:329
    for j = 1:129
        if (~any(chloroL(j,i)))
            tickyboiL = [tickyboiL i];
            break;
        else
            antiTickL = [antiTickL i];
            break;
        end
    end
end

%% Create Arrays of chl-a and time which have no zero or NaN data
chloroNew = chloro(:,[antiTick]);
timeNew = time(:,[antiTick]);
chloroNewL = chloro(:,[antiTickL]);
timeNewL = time(:,[antiTickL]);

%% Look


% MLEpArray = [];
% for i=1:depthMeasurements
%     disp(i);
%     MLEp = statsplot2(chloroNew(i,:));
%     MLEpArray = [MLEpArray MLEp];
% end

%%
depthMeasurements = 129;
clc;
% Eulerian
for i=1:depthMeasurements
    disp(i);
    tmp = chloro(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    [MLEp(:,:,i),KSp(:,i),nll(:,i)] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i=3:depthMeasurements
    disp(i);
    tmp = chloroL(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    [~,KSpL(:,i),nllL(:,i)] = statsplot2(tmp,'noplot');
end

%%
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 15 15]);

ax1 = figure;
subplot(1,2,1)
plot(KSp(1,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(KSp(2,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(KSp(3,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(KSp(4,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(KSp(5,:),linspace(0,-2*depthMeasurements,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian KS-Test');

subplot(1,2,2)
plot(KSpL(1,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(KSpL(2,:),linspace(128,-128,depthMeasurements),'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(KSpL(3,:),linspace(128,-128,depthMeasurements),'Color','red','DisplayName','Weibull');
plot(KSpL(4,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle','--','DisplayName','Gamma');
plot(KSpL(5,:),linspace(128,-128,depthMeasurements),'Color','red','LineStyle',':','DisplayName','Exp');
hold off
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian KS-Test');
exportgraphics(ax1,'figures/ks-EulAndLag_1988_2021.png');