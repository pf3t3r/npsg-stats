clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 28 12]);
set(0,'defaultAxesFontSize',12);

%% Load Test Data
chloro = load("datafiles\chloro.mat"','chloro256_n').chloro256_n;
time = load("datafiles\chloro.mat","time256").time256;

%% Remove NaN or Zeros

tickyboi = [];
antiTick = [];

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

%% Create Arrays of chl-a and time which have no zero or NaN data
chloroNew = chloro(:,[antiTick]);
timeNew = time(:,[antiTick]);

%% Look

R = []; p2 = [];
depthMeasurements = 129;

for i = 1:depthMeasurements
%     disp(i);
    [R(:,i),p2(:,i)] = bbvuong(chloroNew(i,:));
end


%% show

ax1 = figure;
subplot(1,2,1)
plot(R(1,:),linspace(0,-2*depthMeasurements,129),'DisplayName','Normal-Lognormal');
hold on
plot(R(2,:),linspace(0,-2*depthMeasurements,129),'DisplayName','Normal-Weibull');
xline(0,'HandleVisibility','off');
hold off
legend('Location','best');
xlabel('LLR');
ylabel('Depth [m]');
title('Log-Likelihood Ratio');

subplot(1,2,2)
plot(p2(1,:),linspace(0,-2*depthMeasurements,129),'DisplayName','Normal-Lognormal');
hold on
plot(p2(2,:),linspace(0,-2*depthMeasurements,129),'DisplayName','Normal-Weibull');
hold off
xlabel('p-value');
legend('Location','best');
title('Vuong`s Test p-value');

exportgraphics(ax1,'figures/vuong-Eulerian-nl-nw_1988_2021.png');
% [R,p2] = bbvuong(chloroNew(55,:));