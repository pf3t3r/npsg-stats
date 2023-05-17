function [ks, obs, depth2,Sk] = ksOfBinnedCon(X, p, binning, threshold)
%ksOfBinnedCon find the KS statistic
% INPUTS:
% X = substance concentration,
% p = binned pressure [dbar]
% binning = range of depths to bin values (default=10)
% threshold = no. of values needed for us to consider results (default=100)
% OUTPUTS: 
% ks = KS test for five distributions 
% obs = observations per depth
% depth2 = array of depths above threshold (=100)

if nargin <4
    threshold = 100;
end

if nargin <3
    binning = 10;
end

if binning == 5
    ks = nan(5,40); obs = nan(40,1); n = 40; depth = 5:5:200;
elseif binning == 10
    ks = nan(5,20); obs = nan(20,1); n = 20; depth = 5:10:200;
else
    msg = 'Binning input not valid. Must be either 5 or 10 dbar.';
    error(msg);
end

for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = X(p==i);
    % apply KS test to X_i
    % change limit below to >3 to fix error with picoeu -> may change other
    % results
    if length(X_i) > 3
        [~,ks(:,i),~] = statsplot2(X_i,'noplot');
        Sk(i) = skewness(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

for i = 1:n
    if obs(i) < threshold
        ks(:,i) = nan;        
    end
end

tmp = [];
for i = 1:n
    if ~isnan(sum(ks(:,i)))
        tmp = [tmp i];
    end
end
depth2 = depth(tmp);
Sk = Sk(tmp);

ks = ks(:,~all(isnan(ks)));

end