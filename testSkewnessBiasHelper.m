function [mN,mL,sN,sL] = testSkewnessBiasHelper(s,runs,dispFig)
% INPUTS
% s = sample size
% OUTPUTS
% sn = skewness of normally-distributed random data
% sl = skewness of lognormally-distributed random data

if nargin < 3
    dispFig = false;
end

if nargin < 2
    runs = 100;
end

% NORMAL case
sn = nan(runs,1);
for i = 1:runs
    a = randn(s,1);
    sn(i) = skewness(a);
end

% LOGNORMAL case
sl = nan(runs,1);
for i = 1:runs
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

% Get MEAN and STD of runs of skewness
mN = mean(sn);
mL = mean(sl);
sN = std(sn);
sL = std(sl);

end