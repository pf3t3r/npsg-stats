clear;clc;close all;

chla = importdata("data\L0\chlaCtd250dbar.txt").data(:,4);
chla(chla==-9 | chla<=0) = nan;
% N = makedist("Normal","mu",1,"sigma",0.5);
% L = makedist("Lognormal","mu",-2.5,"sigma",0.5);
% G = makedist("Gamma","a",2,"b",0.5);
% W = makedist("Weibull","a",1,"b",3);
% 
% figure;
% plot(N); hold on
% plot(L);
% plot(G);
% plot(W); 
% histogram(chla);
% hold off
% legend("Normal","Lognormal","Gamma","Weibull");

%%
N = fitdist(chla,"Normal");
L = fitdist(chla,"Lognormal");
G = fitdist(chla,"Gamma");
W = fitdist(chla,"Weibull");

figure;
plot(N); hold on
plot(L);
plot(G);
plot(W); hold off
legend("Normal","Lognormal","Gamma","Weibull");
