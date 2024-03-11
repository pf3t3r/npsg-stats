function [ax,ks,obs] = L0_helper(tmp)
%L0_helper
%INPUT: tmp = text-file with pressure, bottle concentration, and bottle ID;
%OUTPUT: ks = K-S p-value, obs = no. of observations, 
% SAVE and OUTPUT bottle ID: reminder!

threshold = 50; tmpT = ""; alphaKs = 0.05;

if isstruct(tmp)
    pIn = tmp.data(:,4);
    X = tmp.data(:,5);
else
    pIn = tmp(:,4);
    X = tmp(:,5);
end
n = length(pIn);

%% Bin
pB = discretize(pIn,0:10:200);
pX = nan(n,1);

for i = 1:n
    pX(i) = pB(i)*10 - 5;
end

%% K-S: bottle chla

n2 = 20;
ks = nan(5,n2);
sk = nan(1,n2); ku = nan(1,n2);
rV = nan(10,n2); pV = nan(10,n2);
obs = nan(1,n2);

% Calculate KS p-value, skewness, kurtosis
for i = 1:n2
    % find concentration X_i at binned pressure i
    X_i = X(pB==i);
    % remove negative or null values
    X_i(X_i <= 0) = [];
    % apply KS test to X_i (only when at least 3 values at binned pressure)
    if length(X_i) > 3
        [~,ks(:,i),~] = statsplot2(X_i,'noplot');
        [rV(:,i),pV(:,i)] = bbvuong(X_i);
        sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

ks(:,obs<threshold) = nan;
sk(obs<threshold) = nan;
ku(obs<threshold) = nan;

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% Negative Distributions
skLognN = -skLogn;
kuLognN = kuLogn;

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

pXX = 5:10:195;
ax = figure;
subplot(1,5,1)
barh(obs);
xline(50);
ylim([0.5 20.5]);
yticks(0.5:2:20.5);
yticklabels(0:20:200);
set(gca,"YDir","reverse");
xlabel("# Observations",FontSize=15);
ylabel("P [dbar]",FontSize=15);

subplot(1,5,[2 3])
plot(ks(2,:),pXX,'+--','Color','#1f78b4',LineWidth=1.5,MarkerSize=5,HandleVisibility='off');
xline(0.05,DisplayName='\alpha'); ylim([0 200]);
set(gca,'YDir','reverse');
yticklabels([]);
grid minor;
legend(FontSize=15);
xlabel('K-S $p$-value',Interpreter='latex',FontSize=15);

clr = 1:1:length(pXX);
subplot(1,5,[4.1 5])
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
hold on
plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
for i = 1:n2
    if ks(2,i) < alphaKs
        plot(sk(i),ku(i),Marker="pentagram",Color='k',HandleVisibility='off',MarkerSize=9);
    else
        plot(sk(i),ku(i),Marker="o",Color=[0.8 0.8 0.8],HandleVisibility='off');
    end
end
scatter(sk,ku,24,clr,"filled","o",HandleVisibility="off");
hold off
colormap(gca,cbrewer2("RdYlBu"));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:1:length(pXX);
cbar.TickLabels = pXX;
cbar.Label.String = "P [dbar]";
legend(Location="best",FontSize=15);
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
ylabel('Kurtosis',FontSize=15); xlabel('Skewness',FontSize=15);

% sgtitle("L0: Bottle Chl-a");
% exportgraphics(ax,"figures/L0/bottle/chla" + tmpT + ".png"); clear ax;

end