function [mN,mL,sN,sL,seN,seL] = testKurtosisBiasHelper(s,runs,dispFig)
% INPUTS
% s = sample size
% OUTPUTS
% kn = kurtosis of normally-distributed random data
% kl = kurtosis of lognormally-distributed random data

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

% Get MEAN and STD of runs of kurtosis
mN = mean(kn);
mL = mean(kl);
sN = std(kn);
sL = std(kl);

% Standard Error
seN = sN/sqrt(length(runs));
seL = sL/sqrt(length(runs));

end