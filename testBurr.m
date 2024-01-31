clear;clc;close all; addpath("baroneRoutines\");
data = importdata("data\L0\chlaCtd50dbar.txt").data;

chla50 = data(:,4);
chla50(chla50==-9) = nan;

pd1 = fitdist(chla50,"Lognormal");
pd2 = fitdist(chla50,"Burr");
pd3 = fitdist(chla50,"Beta");

p_beta = pdf(pd3,sortrows(chla50));
p_burr = pdf(pd2,sortrows(chla50));
p_logn = pdf(pd1,sortrows(chla50));

figure;
yyaxis left
plot(sortrows(chla50),p_burr,'-',sortrows(chla50),p_logn,'-.',sortrows(chla50),p_beta,'.');
hold on
yyaxis right
histogram(chla50);
hold off
legend('Burr','Lognormal','Beta','Data');

%% K-S Comparison

tmp = chla50;
tmp(isnan(tmp)) = [];
[~,ks,~] = statsplot2(tmp,'noplot');

% x-dim of CDF
x_cdf = linspace(min(tmp)-2*std(tmp),max(tmp)+2*std(tmp),2000);

% Burr
paramsBurr = mle(tmp,'distribution','Burr','Alpha',0.32);
y_cdf_burr = cdf('burr',x_cdf,paramsBurr(1),paramsBurr(2),paramsBurr(3));
[hBu,ksBu] = kstest(tmp,[x_cdf' y_cdf_burr']);

% Beta
paramsBeta = mle(tmp,'distribution','Beta','Alpha',0.32);
y_cdf_beta = cdf('Beta',x_cdf,paramsBeta(1),paramsBeta(2));
[hBe,ksBe] = kstest(tmp,[x_cdf' y_cdf_beta']);

%%

sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% BETA

for j = 1:2
    r = randi([1 50],4,1);
    gam = linspace(r(1),r(2),1e5);
    eta = linspace(r(3),r(4),1e5);
    for i = 1:length(gam)
        skBeta(i) = (2*(eta(i) - gam(i))*sqrt(gam(i) + eta(i) + 1) ) ./ ( (gam(i) + eta(i) + 2)*sqrt(gam(i)*eta(i)) );
        kuBeta(i) = 3*(eta(i) + gam(i) + 1) * ( 2*(eta(i) + gam(i))^2 + gam(i)*eta(i)*(gam(i) + eta(i) - 6) )  ./ ...
            ( eta(i)*gam(i)*(eta(i) + gam(i) + 2)*(eta(i) + gam(i) + 3) );
    end
    SKB(j,:) = skBeta;
    KUB(j,:) = kuBeta;
end


figure;
% plot(skLogn,kuLogn);
% hold on
% scatter(skBeta,kuBeta);
scatter(SKB,KUB);
% hold off
% ylim([0 10]); xlim([0 2.5]);