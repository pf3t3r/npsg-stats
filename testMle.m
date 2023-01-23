clear; clc; close all;

load carbig.mat;

%%

figure
histogram(MPG);

%%

% normal
phat = mle(MPG);
des_signal = log(normpdf(MPG,phat(1),phat(2)));


figure
plot(des_signal);
%%
% lognormal
phatln = mle(MPG,'distribution','Lognormal');
%%
% weibull
phatWe = mle(MPG,'distribution','Weibull');
%%
% exponential
phatEx = mle(MPG,'distribution','Exponential');