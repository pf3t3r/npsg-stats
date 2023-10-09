function [ax,p,ks,obs,sk,ku] = L2_helper(tmp,maxMld,dcm,threshold)
%%L2_helper: this function makes the calculation of KS p-values, skewness,
%%and kurtosis a little more efficient for L2 (sub-mixed layer region that
% is centred on the DCM). 
% INPUTS
% tmp: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% dcm: pressure of deep chlorophyll maximum (by cruise)
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
if nargin < 4
    threshold = 50;
end

id = num2str(tmp.data(:,1));
p = tmp.data(:,4);
c = tmp.data(:,5);
clear tmp;

% 2. Extract data beneath ML, centre around DCM
[idSubml,pSubml,cSubml] = extractSMLC(id,p,c,maxMld);

% 3. Calculate KS p-value, skewness, kurtosis
[pr,ks,obs,sk,ku,~,~,~] = ksOfLagrangian(idSubml,pSubml,dcm,cSubml,threshold);

% 4. Plot results
ax = figure;
plotKs2(pr,ks,obs,sk,ku,pr(1),pr(end),threshold);

end