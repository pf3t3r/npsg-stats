clear;clc;close all;

N = makedist("Normal","mu",1,"sigma",0.5);
L = makedist("Lognormal","mu",0,"sigma",0.5);
G = makedist("Gamma","a",2,"b",0.5);
W = makedist("Weibull","a",1,"b",3);

figure;
plot(N); hold on
plot(L);
plot(G);
plot(W); hold off
legend("Normal","Lognormal","Gamma","Weibull");