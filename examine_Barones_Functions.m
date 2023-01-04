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

load("datafiles\chloro.mat");

[testR,testp2] = bbvuong(chloro2D(1:100,1));
% normal dist fits best here....

% [testR_t,testp2_t] = bbvuong(chloro2D(50,:));

% there is an MLE function in MATLAB!
% [phat,pci] = mle(chloro2D(1:100,1));

%% Try: hotFTPextract

cruise = 210;
filename = 'datafiles\ctd_iso' + string(cruise) + '.mat';

if isfile(filename)
    disp('CTD and ISO already created for that cruise');
else
    disp('Extracting...');
    [ctd,iso] = hotFTPextract(cruise);
    save(filename,"ctd","iso");
end

%% Try: RunMedian

testDataRun = RunMedian(chloro2D(:,1),11);

%% plot the above vs original
figure
plot(chloro2D(:,1),-p_grid(:,1));
hold on
plot(testDataRun,-p_grid(:,1));
hold off

%% Try: statsplot2

[MLEp,KSp,nll] = statsplot2(chloro2D(1:100,1));

%%
clear MLEp KSp nll

%% kstest

h_norm = kstest(chloro2D(:,1));
% if testh = 1, then it means that the distribution is normal (?)
