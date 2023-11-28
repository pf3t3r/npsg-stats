function [ax,pL,ks,obs,sk,ku,pV,rV,tr] = L2_helper_FLUORO(X,pIn,maxMld,dcm)
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

threshold = 50;         % Default threshold = 50 [Mishra et al (2019), 
                        % Ghasemi & Zahediasl (2012), Ahad et al (2011)]
n = length(pIn);        % Depth range
nT = length(X(1,:));    % Time range


% 1. Extract data beneath ML

pSubml = nan(n,nT);
xSubml = nan(n,nT);
for i = 1:n-1
    for j = 1:nT
        if pIn(i) >= maxMld(j)
            pSubml(i,j) = pIn(i+1);
            xSubml(i,j) = X(i+1,j);
        end
    end
end


% 2. Convert p to Lagrangian Coordinates

pL = nan(n,nT);
for i = 1:nT
    pL(:,i) = pSubml(:,i) - dcm(i);
end


% 3. Find KS p-value, Vuong p-value and LLR, skewness, and kurtosis for 
% each depth.

l1 = min(min(pL));
l2 = max(max(pL));

range = l1:2:l2;
rangeLen = 1:1:length(range);
n2 = length(range);

ks = nan(5,n2);
rV = nan(10,n2);
pV = nan(10,n2);
obs = nan(1,n2);
sk = nan(1,n2);
ku = nan(1,n2);
%adH = nan(1,n2); adP = nan(1,n2); % optional parameters for A-D Test

for i = rangeLen
    tmp = xSubml(pL==range(i));
    disp(i);
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obs(i) = length(tmp);
    if length(tmp) > 3
        [~,ks(:,i),~,~,~,~] = statsplot2(tmp,'noplot');
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
        %[adH(i),adP(i)] = lillietest(log(tmp),MCTol=1e-2); % optional Lilliefors Test
        sk(i) = skewness(tmp);
        ku(i) = kurtosis(tmp);
    end
end

% Set results = nan if they do not meet the threshold
for i = rangeLen
    if obs(i) < threshold
        ks(:,i) = nan;
        ku(i) = nan;
        pV(:,i) = nan;
        rV(:,i) = nan;
        sk(i) = nan;
    end
end


% 4. Compare all Vuong Test LLR results to give a single best distribution
% for each depth. Comment / uncomment only one of 4.a or 4.b.

vuongRes = zeros(1,n2);
annot = strings(1,n2);
anClr = strings(1,n2);
anClr(cellfun(@isempty,anClr)) = '#FFFFFF';
rV(isnan(rV)) = 0;

% % 4.a. Vuong: Normal vs Lognormal vs Weibull vs Gamma
% for i = 1:n2
%     if rV(1,i) & rV(2,i) & rV(3,i) > 0
%         vuongRes(i) = 1;
%         annot(i) = "Normal";
%         anClr(i) = '#a6cee3';
%     elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
%         vuongRes(i) = 2;
%         annot(i) = "Lognormal";
%         anClr(i) = '#1f78b4';
%     elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
%         vuongRes(i) = 3;
%         %annot(i) = "";
%         annot(i) = "Weibull";
%         anClr(i) = '#b2df8a';
%     elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
%         vuongRes(i) = 4;
%         %annot(i) = "";
%         annot(i) = "Gamma";
%         anClr(i) = '#33a02c';
%     end
% end
rV(rV==0) = nan;

% 4.b. Vuong: Normal Vs. Lognormal Only
for i = 1:n2
    if rV(1,i) > 0 
        vuongRes(i) = 1;
        annot(i) = "Normal";
        anClr(i) = '#a6cee3';
    elseif rV(1,i) < 0
        vuongRes(i) = 2;
        annot(i) = "Lognormal";
        anClr(i) = '#1f78b4';
    else
        annot(i) = "";
    end
end

% 5. Generate theoretical skewness-kurtosis curves for the...
% Lognormal family,
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% Gamma family, and
kTh = linspace(0.2,5000,10000);
for i = 1:length(kTh)
    skGam(i) = 2/sqrt(kTh(i));
    kuGam(i) = 6/kTh(i) + 3;
end

% Weibull family
kWbl = linspace(0,5,10000);
for i = 1:length(kWbl)
    skWbl(i) = ( gamma(1 + 3/kWbl(i)) - 3*gamma(1 + 1/kWbl(i))*gamma(1 + 2/kWbl(i)) + 2*(gamma(1 + 1/kWbl(i)))^3 ) ./ ...
        ( gamma(1 + 2/kWbl(i)) -  (gamma(1 + 1/kWbl(i)))^2 )^(3/2);
    kuWbl(i) = ( gamma(1 + 4/kWbl(i)) - 4*gamma(1 + 1/kWbl(i))*gamma(1 + 3/kWbl(i)) + 6*( (gamma(1 + 1/kWbl(i)) )^2)*gamma(1 + 2/kWbl(i)) - 3*( (gamma(1 + 1/kWbl(i)))^4 ) ) ./ ...
       ( gamma(1 + 2/kWbl(i)) - ( gamma(1 + 1/kWbl(i)) )^2 )^2;
end


% 6. Plot everything.

ax = figure;

subplot(1,6,1)
barh(obs(rangeLen(1):rangeLen(end)),'FaceColor','#a6cee3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylabel('Pressure [dbar]');
set(gca,"YTick",1:5:n2,"YTickLabel",range(1):10:range(end));
ylim([rangeLen(21) rangeLen(end-30)]);
title('No. of Observations');

subplot(1,6,[2 3])
plot(ks(1,:),range,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
hold on
plot(ks(2,:),range,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
% plot(ks(3,:),range,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
% plot(ks(4,:),range,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
hold off
grid minor;
ylim([l1+40 l2-60]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('$p$-value','Interpreter','latex');
title('K-S Test');

zzs = 0.25*ones(n2,1);
subplot(1,6,4)
for i = 1:n2
    text(zzs(i),range(i),annot(i),FontSize=8,Color=anClr(i));
end
hold on
for i = 1:n2
    if pV(1,i) > 0.05
        scatter(pV(1,i),range(i),[],"black");
    end
end
hold off
ylim([l1+40 l2-60]); 
set(gca,'YDir','reverse');
title('Vuong LLR');

tmp = [];
for i = 1:n2
    if ~isnan(sum(ks(:,i)))
        tmp = [tmp i];
    end
end
tr2 = range(tmp);
sk2 = sk(tmp);
ku2 = ku(tmp);
clear tmp;

kurtLimB = 10; skewLimA = 0; skewLimB = 2.5;
if max(ku2) > 10 & min(sk2) < 0
    kurtLimB = max(ku2) + 1;
    skewLimA = min(sk2) - 0.1;
    skewLimB = max(sk2) + 0.1;
elseif max(ku2) > 10
    kurtLimB = max(ku2) + 1;
    skewLimB = max(sk2) + 0.1;
elseif min(sk2) < 0 
    skewLimA = min(sk2) - 0.1;
elseif max(sk2) > 2.5
    kurtLimB = max(ku2) + 1;
    skewLimB = max(sk2) + 0.1;
else 
    kurtLimB = 10;
    skewLimA = 0;
    skewLimB = 2.5;
end

subplot(1,6,[5 6])
clr = 1:1:length(tr2);
scatter(sk2,ku2,24,clr,"filled","o",HandleVisibility="off");
colormap(gca,flipud(colormap("hot")));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:10:length(tr2);
cbar.TickLabels = tr2(1):20:tr2(end);
cbar.Label.String = "P [dbar]";
hold on
plot(skLogn,kuLogn,'DisplayName','Logn.','Color','k',LineWidth=2);
plot(skGam,kuGam,'DisplayName','Gam.','Color','#a6cee3',LineWidth=2);
plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#1f78b4',LineStyle=':',LineWidth=2);
scatter(2,9,'DisplayName','Exp.',Marker='+',LineWidth=1);
scatter(0,9/5,'DisplayName','Uni.',Marker='o',LineWidth=1);
scatter(0,3,'DisplayName','Norm.',Marker='*',LineWidth=1);
scatter(0,21/5,'DisplayName','Logi.',Marker='.',LineWidth=1);
scatter(1.1395,5.4,'DisplayName','LEV',Marker='x',LineWidth=1);
hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness'); ylabel('Kurtosis');
lgd = legend('Location','best');
title(lgd,'Distributions');
title('Skewness vs. Kurtosis');

sk = sk2; ku = ku2; tr = tr2;

end