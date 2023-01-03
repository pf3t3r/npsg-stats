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

% update the check to include cruise number in filename?
if isfile("datafiles\ctd_iso.mat")
    disp('CTD and ISO already created');
else
    disp('Extracting...');
    [ctd,iso] = hotFTPextract(cruise);
    filename = 'datafiles\ctd_iso' + string(cruise);
    save(filename,"ctd","iso");
end