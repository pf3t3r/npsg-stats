% File to open and examine HOTS chloropigment data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',16);

addpath("baroneRoutines\");

%% Try Barone's Functions: bbvuong

load("datafiles\chloro.mat","chloro200",'pgrid200');

[testR,testp2] = bbvuong(chloro200(1:100,1));
% normal dist fits best here....

% [testR_t,testp2_t] = bbvuong(chloro2D(50,:));

% there is an MLE function in MATLAB!
% [phat,pci] = mle(chloro2D(1:100,1));

%% Try: hotFTPextract

cruise = 300;
filename = 'datafiles\ctd_iso' + string(cruise) + '.mat';

if isfile(filename)
    disp('CTD and ISO already created for that cruise');
else
    disp('Extracting...');
    [ctd,iso] = hotFTPextract(cruise);
    save(filename,"ctd","iso");
end

%% isopyncnal view of DCM
% plot sig vs. f

load(filename);

% mean fluorescence for that cruise
f_mean = mean(ctd.f,2,'omitnan');
sig_mean = mean(ctd.sig,2,"omitnan");

figureName = 'figures/isopycnal_DCM_' + string(cruise) + '.png';

ax1 = figure;
scatter(ctd.f,-ctd.sig,'Marker','.','MarkerEdgeColor',[0.6 0.6 0.6],'HandleVisibility','off');
hold on
plot(f_mean,-sig_mean,'Color','red','LineWidth',1.5,'DisplayName','Mean Fluorescence');
hold off
legend();
xlabel('Chloropigment [ug L^{-1}]');
ylabel('Potential Density Anomaly \sigma_0 [kg m^{-3}]');
title(sprintf('Isopycnal View of DCM: \\sigma_0 vs fluorescence (Cruise %s)',string(cruise)));
exportgraphics(ax1,figureName);

%% Try: RunMedian

testDataRun = RunMedian(chloro200(:,1),11);

%% plot the above vs original
figure
plot(chloro200(:,1),-pgrid200(:,1));
hold on
plot(testDataRun,-pgrid200(:,1));
hold off

%% Try: statsplot2

[MLEp,KSp,nll] = statsplot2(chloro200(1:100,1));

%%
clear MLEp KSp nll

%% kstest

h_norm = kstest(chloro200(:,1));
% if testh = 1, then it means that the distribution is normal (?)
