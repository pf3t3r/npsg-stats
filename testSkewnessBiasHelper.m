function [sn,sl] = testSkewnessBiasHelper(s,dispFig)
% INPUTS
% s = sample size
% OUTPUTS
% sn = skewness of normally-distributed random data
% sl = skewness of lognormally-distributed random data

if nargin < 2
    dispFig = false;
end

% NORMAL case
sn = nan(100,1);
for i = 1:100
    a = randn(s,1);
    sn(i) = skewness(a);
end

% LOGNORMAL case
sl = nan(100,1);
for i = 1:100
    b = lognrnd(1,0.3,[s 1]);
    sl(i) = skewness(b);
end

if dispFig == true
    ax = figure;
    plot(sn,DisplayName='Random Normal');
    hold on
    plot(sl,DisplayName='Random Lognormal');
    hold off
    legend();
    title('Skewness of randomly-generated arrays',sprintf('Sample Size = %d',s));
    str = "t" + s + "_.png";
    exportgraphics(ax,'figures/skewBias/' + str); clear ax;
end

end