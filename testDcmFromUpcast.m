clc; close all; clear;
addpath("C:\Users\pfarrell\OneDrive - NIOZ\Documenten\hotFTPupcast");
%% Extract DCM from upcast

for i = 1:10
    [pcm(i),sigcm(i),~] = hotFTPupcast(300,i);
end

%% Compare with hotFTPextract dcm

dcmE = load("dcm.mat").dcm;

%% CRN 300

ax = figure;
plot(dcmE(9270:9279,2),dcmE(9270:9279,3),DisplayName='EXTRACT dcm');
hold on
plot(1:1:10,pcm,DisplayName='UPCAST dcm');
hold off
legend(); xlabel('Cast'); ylabel('Pressure [dbar]');
set(gca,"YDir","reverse");
title('HOT-300: Comparison of DCM Calculation Routines');
exportgraphics(ax,'figures/dcmComp/hot300.png');