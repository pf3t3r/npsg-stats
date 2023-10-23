function [ax,p,ks,obs,Sk,Ku,sd,rV,pV] = L1_helper(tmp,maxMld,threshold)
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

% Set default threshold
% Default threshold of 50 based on findings of Mishra et al (2019), Ghasemi
% & Zahediasl (2012), and Ahad et al (2011).
if nargin < 3
    threshold = 50;
end

pIn = tmp.data(:,4);
cIn = tmp.data(:,5);
idIn = tmp.data(:,1);
clear tmp;

% 2. Extract data in ML
[idOut,pOut,cOut] = extractMldVals(idIn,pIn,cIn,maxMld);

% 3. Bin data
[~,pOutB,cOutB,~,~] = cleanAndBin(pOut,cOut,idOut');

% 4. Calculate KS p-value, skewness, kurtosis
[ks,obs,p,Sk,Ku,sd,rV,pV] = ksOfBinnedCon(cOutB,pOutB,10,threshold);

% 4.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = nan(1,length(p));
for i = 1:length(p)
    if rV(1,i) & rV(2,i) & rV(3,i) > 0
        disp('Normal');
        vuongRes(i) = 1;
    elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
        disp('Lognormal');
        vuongRes(i) = 2;
    elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
        disp('Weibull');
        vuongRes(i) = 3;
    elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
        disp('Gamma');
        vuongRes(i) = 4;
    end
end

% 5. Plot results
ax = figure;
plotKs(p,ks,obs,Sk,Ku,0.5,20.5,true,threshold,vuongRes,rV,pV);

end