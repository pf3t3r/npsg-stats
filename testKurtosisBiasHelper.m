function [mN,mL,pN,pL] = testKurtosisBiasHelper(s,runs,dispFig)
% INPUTS
% s = sample size
% runs = no. of random arrays generated in for loop
% dispFig = set true if you want to see the figure, otherwise false
% OUTPUTS
% kn = kurtosis of normally-distributed random data
% kl = kurtosis of lognormally-distributed random data
% mN = mean kurtosis of normally-distributed random data
% mL = mean kurtosis of lognormally-distributed random data
% pN = lower and upper percentile of the normal kurtoses
% pL = lower and upper percentile of the lognormal kurtoses


if nargin < 3
    dispFig = false;
end

if nargin < 2
    runs = 100;
end

% NORMAL case
kn = nan(runs,1);
for i = 1:runs
    a = randn(s,1);
    kn(i) = kurtosis(a);
end

% LOGNORMAL case
kl = nan(runs,1);
for i = 1:runs
    b = lognrnd(1,0.3,[s 1]);
    kl(i) = kurtosis(b);
end

if dispFig == true
    ax = figure;
    plot(kn,DisplayName='Random Normal');
    hold on
    plot(kl,DisplayName='Random Lognormal');
    hold off
    legend();
    title('Kurtosis of randomly-generated arrays',sprintf('Sample Size = %d',s));
    str = "test" + s + "bias.png";
    exportgraphics(ax,'figures/kurtBias/' + str); clear ax;
end

% Get MEAN and PRCTILE of runs of kurtosis
mN = mean(kn);
mL = mean(kl);
% pN = std(kn);
% pL = std(kl);
pN = [prctile(kn,16) prctile(kn,84)];
pL = [prctile(kl,16) prctile(kl,84)];

% % Standard Error
% seN = pN/sqrt(length(runs));
% seL = pL/sqrt(length(runs));

end