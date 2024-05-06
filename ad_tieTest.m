clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot,'defaultFigureUnits','centimeters','defaultFigurePosition',[3 3 28 15]);

logAxes = true;

tmp = importdata('data/L0/hplcChla_88-21_200.txt');
[ax,ks,~,pB,X] = L0_helper(tmp,50,"ad",logAxes);
% sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");

testX = X(pB==1);

%%
randNum = 0.00099*randi([-1000 1000],length(testX),1);
testX2 = testX + randNum;

figure;
subplot(2,1,1)
histogram(testX);
hold on
subplot(2,1,2)
histogram(testX2);
hold off

[adh,adp] = adtest(testX,"Distribution","logn");
[adh2,adp2] = adtest(testX2,"Distribution","logn");