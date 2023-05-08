function [ks, obs, depth2] = ksOfBinnedCon(X, p, binning)
%funshit2 find the KS statistc
% INPUTS:
% X = substance concentration,
% p = binned pressure [dbar]
% n = length of data set
% OUTPUTS: 
% ks = KS test for five distributions 
% obs = observations per depth
% depth2 = array of depths above threshold (=100)

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
    % apply KS test to chla_i
    if length(X_i) > 1
        [~,ks(:,i),~] = statsplot2(X_i,'noplot');
    end
    obs(i) = length(X_i);
    clear X_i;
end

for i = 1:n
    if obs(i) < 100
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

ks = ks(:,~all(isnan(ks)));

end