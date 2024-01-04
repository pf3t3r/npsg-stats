clear; clc; close all;

x = ncinfo("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc");

t = ncread("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc","TIME");
T = ncread("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc","TEMP");
Sp = ncread("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc","PSAL");
p = ncread("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc","PRES");

%%
figure;
subplot(1,3,1)
scatter(Sp,p,[],[0.4 0.4 0.4],'.');
subplot(1,3,2)
scatter(T,p,[],[0.4 0.4 0.4],'.');
subplot(1,3,3)
scatter(Sp,T,[],[0.4 0.4 0.4],'.');

%%
figure
subplot(3,1,1)
scatter(t,Sp,[],[0.4 0.4 0.4],'.');
subplot(3,1,2)
scatter(t,T,[],[0.4 0.4 0.4],'.');


%%
sam = Sp;
% sam = randsample(Sp,300);

% [h,p] = adtest(sam,"Distribution","logn");
% adtest doesn't work for large sample sizes

[h,p] = lillietest(sam,"Distr","norm");

sk = skewness(Sp);
ku = kurtosis(Sp);