% Generate and plot the theoretical curves of skewness and kurtosis for a
% selection of empirical distributions.

clc; close all; clear;

%% Generate theoretical skewness and kurtosis curves.

% Lognormal family
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% Gamma family
kTh = linspace(0.04,3000,1500000);
for i = 1:length(kTh)
    skGam(i) = 2/sqrt(kTh(i));
    kuGam(i) = 6/kTh(i) + 3;
end

% Weibull family
kWbl = linspace(0.1,3.5,10000);
for i = 1:length(kWbl)
    skWbl(i) = ( gamma(1 + 3/kWbl(i)) - 3*gamma(1 + 1/kWbl(i))*gamma(1 + 2/kWbl(i)) + 2*(gamma(1 + 1/kWbl(i)))^3 ) ./ ...
        ( gamma(1 + 2/kWbl(i)) -  (gamma(1 + 1/kWbl(i)))^2 )^(3/2);
    kuWbl(i) = ( gamma(1 + 4/kWbl(i)) - 4*gamma(1 + 1/kWbl(i))*gamma(1 + 3/kWbl(i)) + 6*( (gamma(1 + 1/kWbl(i)) )^2)*gamma(1 + 2/kWbl(i)) - 3*( (gamma(1 + 1/kWbl(i)))^4 ) ) ./ ...
       ( gamma(1 + 2/kWbl(i)) - ( gamma(1 + 1/kWbl(i)) )^2 )^2;
end

% Add negative curves
skLognN = -skLogn;
kuLognN = kuLogn;
skGamN= -skGam;
kuGamN = kuGam;
skWblN = -skWbl;
kuWblN = kuWbl;

%% Plot skewness-kurtosis curves.

figure;
scatter(0,3,[],[0.6509803921568628 0.807843137254902 0.8901960784313725],'DisplayName','Normal',Marker='o',LineWidth=3);
hold on
plot(skLogn,kuLogn,'Color','#1f78b4',LineStyle='--',DisplayName='Lognormal');
plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',HandleVisibility='off');
plot(skWbl,kuWbl,'Color','#b2df8a',LineStyle='-',DisplayName='Weibull');
plot(skWblN,kuWblN,'Color','#b2df8a',LineStyle='-',HandleVisibility='off');
plot(skGam,kuGam,'Color','#33a02c',LineStyle='--',DisplayName='Gamma');
plot(skGamN,kuGamN,'Color','#33a02c',LineStyle='--',HandleVisibility='off');
hold off
ylim([0 100]);
legend();
% xlim ([-2.5 2.5]);