function [kn,kl] = testKurtosisBiasHelper(s,dispFig)
% INPUTS
% s = sample size
% OUTPUTS
% kn = kurtosis of normally-distributed random data
% kl = kurtosis of lognormally-distributed random data

if nargin < 2
    dispFig = false;
end

% NORMAL case
kn = nan(100,1);
for i = 1:100
    a = randn(s,1);
    kn(i) = kurtosis(a);
end

% LOGNORMAL case
kl = nan(100,1);
for i = 1:100
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

end