function [ax,p,pKs,pLil,pAd,obs,Sk,Ku] = L1_compareHelper(tmp,maxMld,threshold)
%%L1_helper_V2: this function allows for the calculation of the p-values
% according to the Lilliefors-corrected K-S test (Lil), the basic K-S test
% (K-S), and the Anderson-Darling test (A-D).
% INPUTS
% tmp: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% OUTPUTS
% ax = figure identifier, needed for labelling and export,
% p = pressures where sufficient measurements exist in mixed layer,
% pKs = KS test p-values at those pressures,
% obs = observations per depth,
% Sk = skewness at depths where pKs is taken,
% Ku = kurtosis at the same depths.

% Default threshold = 50. This is based on recommendations from Mishra et 
% al (2019), Ghasemi & Zahediasl (2012), and Ahad et al (2011).
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
[pKs,pLil,pAd,obs,p,Sk,Ku] = pValueOfBinned(cOutB,pOutB,10,threshold);

% 5. Plot results
ax = figure;
subplot(1,4,1)
barh(obs,'FaceColor','#a6cee3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([0.5 length(obs)+0.5]);
ylabel('Pressure [dbar]');
set(gca,"YTick",1:1:length(5:10:205),"YTickLabel",5:10:205);
title('No. of Observations');

subplot(1,4,2)
plot(pKs(1,:),p,'o-','Color','#a6cee3','DisplayName','K-S','LineWidth',1.5,'MarkerSize',5);
hold on
plot(pLil(1,:),p,'o-','Color','#1f78b4','DisplayName','Lillie','LineWidth',1.5,'MarkerSize',5);
plot(pAd(1,:),p,'o-','Color','#b2df8a','DisplayName','A-D','LineWidth',1.5,'MarkerSize',5);
grid minor;
ylim([0 200]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
title('Normal');

subplot(1,4,3)
plot(pKs(2,:),p,'o-','Color','#a6cee3','DisplayName','K-S','LineWidth',1.5,'MarkerSize',5);
hold on
plot(pLil(2,:),p,'o-','Color','#1f78b4','DisplayName','Lillie','LineWidth',1.5,'MarkerSize',5);
plot(pAd(2,:),p,'o-','Color','#b2df8a','DisplayName','A-D','LineWidth',1.5,'MarkerSize',5);
grid minor;
ylim([0 200]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
title('Lognormal');

subplot(1,4,4)
plot(pKs(3,:),p,'o-','Color','#a6cee3','DisplayName','K-S','LineWidth',1.5,'MarkerSize',5);
hold on
plot(pLil(3,:),p,'o-','Color','#1f78b4','DisplayName','Lillie','LineWidth',1.5,'MarkerSize',5);
plot(pAd(3,:),p,'o-','Color','#b2df8a','DisplayName','A-D','LineWidth',1.5,'MarkerSize',5);
grid minor;
ylim([0 200]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
title('Weibull');

end