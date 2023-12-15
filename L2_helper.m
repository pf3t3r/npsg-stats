function [ax,p,ks,obs,sk,ku,sd,rV,pSubml,pV] = L2_helper(tmp,maxMld,dcm,tmpLts,threshold)
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
if nargin < 5
    threshold = 50;
end

if nargin < 4
    tmpLts = [2 22];
end

id = num2str(tmp.data(:,1));
p = tmp.data(:,4);
c = tmp.data(:,5);
clear tmp;

% 2. Extract data beneath ML
[idSubml,pSubml,cSubml] = extractSMLC(id,p,c,maxMld);

% 3. Calculate KS p-value, skewness, kurtosis
% ..., centre around DCM (?)
[pr,ks,obs,sk,ku,sd,rV,pV] = ksOfLagrangian(idSubml,pSubml,dcm,cSubml,threshold);

% 3.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = zeros(1,length(pr));
rV(isnan(rV)) = 0;

% for i = 1:length(pr)
%     %disp(i);
%     if rV(1,i) & rV(2,i) & rV(3,i) > 0
%         %disp('Normal');
%         vuongRes(i) = 1;
%     elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
%         %disp('Lognormal');
%         vuongRes(i) = 2;
%     elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
%         %disp('Weibull');
%         vuongRes(i) = 3;
%     elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
%         %disp('Gamma');
%         vuongRes(i) = 4;
%     end
% end
% rV(rV==0) = nan;

for i = 1:length(pr)
    if rV(1,i)  > 0
        vuongRes(i) = 1;
    elseif rV(1,i) < 0
        vuongRes(i) = 2;
    end
end
rV(rV==0) = nan;

limits = [pr(tmpLts(1)) pr(tmpLts(2))];
obsId = [tmpLts(1) tmpLts(2)];

% 4. Plot results
ax = figure;
plotKs2(pr,ks,obs,sk,ku,limits(1),limits(end),threshold,vuongRes,obsId,pV);

end