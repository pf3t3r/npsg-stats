function [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,threshold)
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
if nargin < 3
    threshold = 100;
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
[ks,obs,p,Sk,Ku] = ksOfBinnedCon(cOutB,pOutB,10,threshold);

% 5. Plot results
ax = figure;
plotKs(p,ks,obs,Sk,Ku,0.5,20.5,true,threshold);

end