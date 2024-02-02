%% Histograms
close all;clc;clear;
tmp = importdata("data\L0\dvchla_94-21_200.txt").data;
p = tmp(:,4); 
X = tmp(:,5);
n = length(p);

pB = discretize(p,0:10:200);
pX = nan(n,1);
for i = 1:n
    pX(i) = pB(i)*10 - 5;
end
clear i n p pB tmp;
%%
NAME = "dvchla-25";
tmpX = X(pX==25);
tmpX(tmpX<=0) = [];

figure;
histfit(tmpX,40,"lognormal");
title(NAME);
legend("data","lognormal");
exportgraphics(gca,"figures/L0/bot/hist/" + NAME + ".png")

%%
[h,p] = adtest(tmpX,"Distribution","logn");
[hl,pl] = lillietest(log(tmpX),"Distr","norm");
% [hc,pc] = jbtest(tmpX);