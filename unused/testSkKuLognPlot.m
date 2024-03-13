clc; close all; clear;

% Script to output the theoretical relationship between skewness and
% kurtosis for the lognormal family of distributions.

sigma = linspace(0,1,1000);

for i = 1:length(sigma)
    sk(i) = (exp(sigma(i)^2) + 2)*(sqrt(exp(sigma(i)^2) - 1));
    ku(i) = exp(4*sigma(i)^2) + 2*exp(3*sigma(i)^2) + 3*exp(2*sigma(i)^2) - 3;
end

% Weibull family: generate theoretical skewness and kurtosis
kWbl = linspace(0.5,100,1000);
gF = 0.2;
for i = 1:length(kWbl)
    skWbl(i) = (gF*(1+3/kWbl(i)) - 3*gF*(1+1/kWbl(i))*gF*(1+2/kWbl(i)) + 2*gF^3*(1+1/kWbl(i))) ./ ...
        (gF*(1+2/kWbl(i)) - gF^2*(1+1/kWbl(i)))^(3/2);
    kuWbl(i) = ( gF*(1+4/kWbl(i)) - 4*gF*(1+1/kWbl(i))*gF*(1+3/kWbl(i)) + 6*gF^2*(1+1/kWbl(i))*gF*(1+2/kWbl(i)) - 3*gF^4*(1+1/kWbl(i))) ./ ...
       ((gF*(1+2/kWbl(i)) - gF^2*(1+1/kWbl(i)))^2);
end

%%
figure;
scatter(skWbl,kuWbl,'.');
% xlim([0 3]);
xlabel('Skewness');
% ylim([0 10]);
ylabel('Kurtosis');