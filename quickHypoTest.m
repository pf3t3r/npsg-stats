close all;clc;clear;
%% Data import / prep.

% Possible File Formats
% y = ncinfo("data\whots\OS_WHOTS_2019_D_M.nc");
uwnd = ncread("data\whots\OS_WHOTS_2019_D_M.nc","UWND");
vwnd = ncread("data\whots\OS_WHOTS_2019_D_M.nc","VWND");
% tmp = importdata('data/L0/hplcChla_88-21_200.txt');

U = sqrt(uwnd.^2 + vwnd.^2);

%% Sample / or not
UU = randsample(U,400);

%% Hypothesis Tests: K-S, Lilliefors, A-D, S-W
% K-S: Most distributions. Depends on CDF function.
% Lilliefors: Normal, Lognormal, Weibull, Exponential, Extreme Value.
% A-D: Normal, Lognormal, Weibull, Exponential, Extreme Value
% S-W: Normal (only).
DIST = "Gamma";

phat = mle(UU,distribution=DIST);
x_cdf = linspace(min(UU)-2*std(UU),max(UU)+2*std(UU),2000);
y_cdf_norm = cdf(DIST,x_cdf,phat(1),phat(2));
[hK,pK] = kstest(UU,[x_cdf' y_cdf_norm']);

% DIST = "norm";
% [hL,pL] = lillietest(log(UU),Distr=DIST);
% [hT,pT] = adtest(UU,"Distribution",DIST);
% [hS,pS] = swtest(UU);

figure;
histfit(UU,[],"weibull");