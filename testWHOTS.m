clear; clc; close all;

% import whots-16 mooring and met station data
x = ncinfo("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc");
y = ncinfo("data\whots\OS_WHOTS_2019_D_M.nc");

t = ncread("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc","TIME");
T = ncread("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc","TEMP");
Sp = ncread("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc","PSAL");
p = ncread("data\whots\OS_WHOTS_201910_D_MICROCAT-120m.nc","PRES");

uwnd = ncread("data\whots\OS_WHOTS_2019_D_M.nc","UWND");
vwnd = ncread("data\whots\OS_WHOTS_2019_D_M.nc","VWND");
tMet = ncread("data\whots\OS_WHOTS_2019_D_M.nc","TIME");
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
figure
subplot(2,1,1)
plot(tMet,uwnd);
xlabel("time"); ylabel("eastward wind (m/s)");
subplot(2,1,2)
plot(tMet,vwnd);
xlabel("time"); ylabel("northward wind (m/s)");

%%
sam = Sp;
% sam = randsample(Sp,300);

% [h,p] = adtest(sam,"Distribution","logn");
% adtest doesn't work for large sample sizes

[h,p] = lillietest(sam,"Distr","norm");

sk = skewness(Sp);
ku = kurtosis(Sp);

%%
% .....[..1.....2.....3....4......5.....6.....7.....8.....9.....10.]
% p2 = [p2_nl;p2_nw;p2_ng;p2_ne;p2_lw;p2_lg;p2_le;p2_wg;p2_we;p2_ge];

[R,p2] = bbvuong(abs(vwnd));

[Ru,p2u] = bbvuong(abs(uwnd));