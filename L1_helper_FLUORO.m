function [ax,tr,ks,obs,sk,ku,sd,rV,pV] = L1_helper_FLUORO(X,pIn,maxMld,threshold)
%%L1_helper: this function makes the calculation of KS p-values, skewness,
%%and kurtosis a little more efficient for L1 (the mixed layer). 
% INPUTS
% tmp: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% OUTPUTS
% ax = figure identifier, needed for labelling and export,
% p = pressures where sufficient measurements exist in mixed layer,
% ks = KS test p-values at those pressures,
% obs = observations per depth,
% Sk = skewness at depths where ks is taken,
% Ku = kurtosis at the same depths.

% Default threshold = 50 based on findings of Mishra et al (2019), Ghasemi
% & Zahediasl (2012), and Ahad et al (2011).
if nargin < 4
    threshold = 50;
end

% for i = 1:length(epN(1,:))
%     test = maxMld(i);
%     disp(test);
% end

n = length(pIn);

copyX = nan(size(X));
for k = 1:length(X(1,:))
    for j = 1:n
        if pIn(j) < maxMld(k)
            copyX(j,k) = X(j,k);
        end
    end
end


% init
ks = nan(5,n); rV = nan(10,n); pV = nan(10,n); sd = nan(n,2); sk = nan(1,n);
ku = nan(1,n); obs = nan(1,n);

for i = 1:n
    tmp = copyX(i,:);
    tmp(isnan(tmp) | tmp<0) = 0;
    tmp(tmp==0) = [];
    if length(tmp) > 3
        [~,ks(:,i),~,~,sd(i,:),~] = statsplot2(tmp,'noplot');
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
        sk(i) = skewness(tmp);
        ku(i) = kurtosis(tmp);
        obs(i) = length(tmp);
    end
end

% tr = pIn;

for i = 1:n
    if obs(i) < threshold
        ks(:,i) = nan;   
        sk(i) = nan;
        ku(i) = nan;
        sd(i,:) = nan;
        rV(:,i) = nan;
        pV(:,i) = nan;
    end
end

tmp = [];
for i = 1:n
    if ~isnan(sum(ks(:,i)))
        tmp = [tmp i];
    end
end
tr = pIn(tmp);
sk = sk(tmp);
ku = ku(tmp);
sd = sd(tmp,:);
rV = rV(:,tmp);
pV = pV(:,tmp);

ks = ks(:,~all(isnan(ks)));

% 4.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = nan(1,length(tr));

% % 4.a.i. Default Case.
% for i = 1:length(tr)
%     if rV(1,i) & rV(2,i) & rV(3,i) > 0
%         disp('Normal');
%         vuongRes(i) = 1;
%     elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
%         disp('Lognormal');
%         vuongRes(i) = 2;
%     elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
%         disp('Weibull');
%         vuongRes(i) = 3;
%     elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
%         disp('Gamma');
%         vuongRes(i) = 4;
%     end
% end

% 4.a.ii. Normal vs Lognormal Case.
for i = 1:length(tr)
    if rV(1,i) > 0
        disp('Normal');
        vuongRes(i) = 1;
    elseif rV(1,i) < 0
        disp('Lognormal');
        vuongRes(i) = 2;
    end
end

% 5. Plot results
ax = figure;
plotKs(tr,ks,obs,sk,ku,0,94,true,threshold,vuongRes,pV,[0 94],true);
% 
% disp(vuongRes);
end