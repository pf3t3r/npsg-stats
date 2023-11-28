function [] = plotKs(tr,ks,obs,sk,ku,obsLimA,obsLimB,EulLan,threshold,vuongRes,pV,limitOveride,fluoOveride)
%plotKs
% INPUT: 
% OUTPUT: 

if nargin < 13
    fluoOveride = false;
end

if nargin < 10
    vuongRes = 0;
end

if nargin < 9
    threshold = 50;
end

if nargin < 8
    EulLan = true;
end

if nargin < 6
    obsLimA = 0.5;
    %obsLimB = length(tr);
    obsLimB = length(obs) + 0.5;
end

if EulLan
    limits = [0 200];
    if obsLimB-obsLimA > 20
        ytix = 5:5:205;
    else
        ytix = 5:10:205;
    end
else
    limits = [-150 100];
    ytix = tr;
end

if nargin >= 12
    limits = limitOveride;
end

if fluoOveride
    ytix = obsLimA:2:obsLimB;
end

n = length(tr);

% Create Annotations for Vuong's Test Results
annot = strings(1,n);
anClr = strings(1,n);
anClr(cellfun(@isempty,anClr)) = '#FFFFFF';

% % Default Case
% for i = 1:n
%     if vuongRes(i) == 1
%         annot(i) = "Normal";
%         anClr(i) = '#a6cee3';
%     elseif vuongRes(i) == 2
%         annot(i) = "Lognormal";
%         anClr(i) = '#1f78b4';
%     elseif vuongRes(i) == 3
%         annot(i) = "Weibull";
%         anClr(i) = '#b2df8a';
%     elseif vuongRes(i) == 4
%         annot(i) = "Gamma";
%         anClr(i) = '#33a02c';
%     elseif vuongRes(i) == 0
%         annot(i) = "";
%     end
% end

% Normal-Lognormal Case
for i = 1:n
    if vuongRes(i) == 1
        annot(i) = "Normal";
        anClr(i) = '#a6cee3';
    elseif vuongRes(i) == 2
        annot(i) = "Lognormal";
        anClr(i) = '#1f78b4';
    elseif vuongRes(i) == 0
        annot(i) = "";
    end
end

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% Gamma family: generate theoretical skewness and kurtosis
kTh = linspace(0.2,5000,10000);
for i = 1:length(kTh)
    skGam(i) = 2/sqrt(kTh(i));
    kuGam(i) = 6/kTh(i) + 3;
end

% Weibull family: generate theoretical skewness and kurtosis
kWbl = linspace(0,5,10000);
for i = 1:length(kWbl)
    skWbl(i) = ( gamma(1 + 3/kWbl(i)) - 3*gamma(1 + 1/kWbl(i))*gamma(1 + 2/kWbl(i)) + 2*(gamma(1 + 1/kWbl(i)))^3 ) ./ ...
        ( gamma(1 + 2/kWbl(i)) -  (gamma(1 + 1/kWbl(i)))^2 )^(3/2);
    kuWbl(i) = ( gamma(1 + 4/kWbl(i)) - 4*gamma(1 + 1/kWbl(i))*gamma(1 + 3/kWbl(i)) + 6*( (gamma(1 + 1/kWbl(i)) )^2)*gamma(1 + 2/kWbl(i)) - 3*( (gamma(1 + 1/kWbl(i)))^4 ) ) ./ ...
       ( gamma(1 + 2/kWbl(i)) - ( gamma(1 + 1/kWbl(i)) )^2 )^2;
end

% Beta family: generate theoretical skewness and kurtosis
muB = linspace(0.5,5,100);
nuB = linspace(0.5,5,100);
alphaB = muB.*nuB; betaB = (1 - muB).*nuB;
for i = 1:length(betaB)
    skBet(i) = 2*(betaB(i) - alphaB(i)) * sqrt(alphaB(i) + betaB(i) + 1) ./ ...
        sqrt(betaB(i)*alphaB(i)) * (betaB(i) + alphaB(i) + 2);
    kuBet(i) = 3*(betaB(i) + alphaB(i) + 1) * ( 2*(betaB(i) + alphaB(i))^2 + alphaB(i)*betaB(i)*(betaB(i) + alphaB(i) - 6) ) ./ ...
        betaB(i)*alphaB(i)*(betaB(i) + alphaB(i) + 2)*(betaB(i) + alphaB(i) + 3);
end

% Loglogistic family: generate theoretical skewness and kurtosis
% Nothings shows up: Maybe I need to tune the shape parameter a bit (??)
% cLgi = linspace(2,4,10000);
% for i = 1:length(cLgi)
%     skLgi(i) = ( 2*(pi^2)*(csc(pi/cLgi(i)))^3 - 6*cLgi(i)*pi*csc(pi/cLgi(i))*csc(2*pi/cLgi(i)) + 3*(cLgi(i)^2)*csc(3*pi/cLgi(i)) ) ./ ...
%         sqrt(pi)*((-pi*(csc(pi/cLgi(i))^2)) + 2*cLgi(i)*csc(2*pi/cLgi(i)))^(3/2);
%     kuLgi(i) = ( -3*pi^2*(csc(pi/cLgi(i))^4) - 12*pi*cLgi(i)^2*csc(pi/cLgi(i))*csc(3*pi/cLgi(i)) + 4*cLgi(i)^3*csc(4*pi/cLgi(i)) + 6*cLgi(i)*pi^2*(csc(pi/cLgi(i))^3)*sec(pi/cLgi(i)) ) ./ ...
%         pi*( -pi*(csc(pi/cLgi(i))^2) + 2*cLgi(i)*csc(2*pi/cLgi(i)) )^2;
% end

subplot(1,6,1)
barh(obs,'FaceColor','#a6cee3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
if fluoOveride
    ylim([obsLimA+1 obsLimB/2 + 1]);
else
    ylim([obsLimA obsLimB]);
end
ylabel('Pressure [dbar]');
set(gca,"YTick",1:1:length(ytix),"YTickLabel",ytix);
title('No. of Observations');

subplot(1,6,[2 3])
plot(ks(1,:),tr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
hold on
plot(ks(2,:),tr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
% plot(ks(3,:),tr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
% plot(ks(4,:),tr,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
hold off
grid minor;
ylim(limits);
set(gca,'YDir','reverse');
legend(Location="best");
xlabel('p-value');
title('K-S p-values');

subplot(1,6,4)
zzs = 0.25*ones(n,1);
for i = 1:n
    %text(zzs(i),tr(i),annot(i),FontSize=8,Color=anClr(i));
    if pV(1,i) > 0.05
        text(zzs(i),tr(i),annot(i),FontSize=8,Color=anClr(i),FontWeight="bold");
    else 
        text(zzs(i),tr(i),annot(i),FontSize=8,Color=anClr(i));
    end
end
% % For Normal-Lognormal Comparison ONLY
% hold on
% pV(1,obs<threshold) = nan;
% for i = 1:n
%     if pV(1,i) > 0.05
%         scatter(pV(1,i),tr(i),[],"black");
%     end
% end
% hold off
grid minor;
set(gca,"Xtick",[]);
ylim(limits); set(gca,'YDir','reverse');
% legend(Location="south");
title('Vuong LLR');

kurtLimB = 10; skewLimA = 0; skewLimB = 2.5;
if max(ku) > 10 & min(sk) < 0
    kurtLimB = max(ku) + 1;
    skewLimA = min(sk) - 0.1;
    skewLimB = max(sk) + 0.1;
elseif max(ku) > 10
    kurtLimB = max(ku) + 1;
    skewLimB = max(sk) + 0.1;
elseif min(sk) < 0 
    skewLimA = min(sk) - 0.1;
elseif max(sk) > 2.5
    kurtLimB = max(ku) + 1;
    skewLimB = max(sk) + 0.1;
else 
    kurtLimB = 10;
    skewLimA = 0;
    skewLimB = 2.5;
end

subplot(1,6,[5 6])
clr = 1:1:length(tr);
scatter(sk,ku,24,clr,"filled","o",HandleVisibility="off");
colormap(gca,flipud(colormap("hot")));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:5:length(tr);
cbar.TickLabels = tr(1):10:tr(end);
cbar.Label.String = "P [dbar]";
hold on
plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#a6cee3',LineStyle='-',LineWidth=1);
plot(skGam,kuGam,'DisplayName','Gam.','Color','#1f78b4',LineStyle='--',LineWidth=1);
plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#b2df8a',LineStyle=':',LineWidth=1);
scatter(2,9,'DisplayName','Exp.',Marker='+',LineWidth=1);
scatter(0,9/5,'DisplayName','Uni.',Marker='o',LineWidth=1);
scatter(0,3,'DisplayName','Norm.',Marker='*',LineWidth=1);
scatter(0,21/5,'DisplayName','Logi.',Marker='.',LineWidth=1);
scatter(1.1395,5.4,'DisplayName','LEV',Marker='x',LineWidth=1);
hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness'); ylabel('Kurtosis');
lgd = legend('Location','best','FontSize',6);
title(lgd,'P [dbar]');
title('Skewness vs. Kurtosis');

end