function [ax,ks,obs] = L0_helper(tmp,threshold,hypTest,logAxis,season)
%L0_helper
%INPUT: tmp = text-file with pressure, bottle concentration, and bottle ID;
%OUTPUT: ks = K-S p-value, obs = no. of observations, 
% SAVE and OUTPUT bottle ID: reminder!

if nargin < 5
    season = 0;
    % 0 = no seasonal analysis, 1 = winter, 2 = spring, 3 = summer,
    % 4 = autumn
end

if nargin < 4
    logAxis = true;
    % true => output p-values in log x-axis, otherwise no log plot.
end

if nargin < 3
    hypTest = 'ks';
end

if nargin <2
    threshold = 50;
end
tmpT = ""; alphaKs = 0.05;

if isstruct(tmp)
    botId = tmp.data(:,2);
    pIn = tmp.data(:,4);
    X = tmp.data(:,5);
else
    botId = tmp(:,2);
    pIn = tmp(:,4);
    X = tmp(:,5);
end
n = length(pIn);
nB = length(botId);
n3 = length(X);

% CONTINUE WORKING ON
if season ~= 0
    botId(botId==-9) = nan;
    botId2 = num2str(botId);
    botMth = nan(n,1);
    for i = 1:n
        tmpX = str2num(botId2(i,1:end-4));
        if ~isnan(tmpX)
            botMth(i) = tmpX;
        end
    end
    winter = nan(n,1); spring = nan(n,1); summer = nan(n,1); autumn = nan(n,1);
    for i = 1:n
        tmpY = botMth(i);
        if (tmpY == 1) || (tmpY == 2) || (tmpY == 3)
            winter(i) = 1;
        end
        if (tmpY == 4) || (tmpY == 5) || (tmpY == 6)
            spring(i) = 1;
        end
        if (tmpY == 7) || (tmpY == 8) || (tmpY == 9)
            summer(i) = 1;
        end
        if (tmpY == 10) || (tmpY == 11) || (tmpY == 12)
            autumn(i) = 1;
        end
    end

    winIds = []; sprIds = []; sumIds = []; autIds = []; 
    for i = 1:n
        if winter(i) == 1
            winIds = [winIds i];
        end
        if spring(i) == 1
            sprIds = [sprIds i];
        end
        if summer(i) == 1
            sumIds = [sumIds i];
        end
        if autumn(i) == 1
            autIds = [autIds i];
        end
    end

    if season == 1
        X = X(winIds);
        pIn = pIn(winIds);
        nSea = length(winIds);
    elseif season == 2
        X = X(sprIds);
        pIn = pIn(sprIds);
        nSea = length(sprIds);
    elseif season == 3
        X = X(sumIds);
        pIn = pIn(sumIds);
        nSea = length(sumIds);
    elseif season == 4
        X = X(autIds);
        pIn = pIn(autIds);
        nSea = length(autIds);
    end
end

%% Bin
pB = discretize(pIn,0:10:200);
pX = nan(n,1);

if season == 0
    for i = 1:n
        pX(i) = pB(i)*10 - 5;
    end
else
    for i = 1:nSea
        pX(i) = pB(i)*10 - 5;
    end
end

%% K-S: bottle chla

n2 = 20;
ks = nan(5,n2);
sk = nan(1,n2); ku = nan(1,n2);
rV = nan(10,n2); pV = nan(10,n2);
obs = nan(1,n2);
ad = nan(1,n2);

% Calculate KS p-value, skewness, kurtosis
for i = 1:n2
    % find concentration X_i at binned pressure i
    X_i = X(pB==i);
    % remove negative or null values
    X_i(X_i <= 0) = [];
    % apply KS test to X_i (only when at least 3 values at binned pressure)
    if length(X_i) > 3
        if strcmp(hypTest,'ks')
            [~,ks(:,i),~] = statsplot2(X_i,'noplot');
        else
            [~,ad(i)] = adtest(X_i,'Distribution','logn');
        end
        [rV(:,i),pV(:,i)] = bbvuong(X_i);
        sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

ks(:,obs<threshold) = nan;
ad(obs<threshold) = nan;
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
xline(threshold);
ylim([0.5 20.5]);
yticks(0.5:2:20.5);
yticklabels(0:20:200);
set(gca,"YDir","reverse");
xlabel("# Observations",FontSize=15);
ylabel("P [dbar]",FontSize=15);

subplot(1,5,[2 3])
xline(0.05,DisplayName='\alpha'); 
hold on
if strcmp(hypTest,'ks')
    plot(ks(2,:),pXX,'+--','Color','#1f78b4',LineWidth=1.5,MarkerSize=5,HandleVisibility='off');
    xlabel('K-S $p$-value',Interpreter='latex',FontSize=15);
else
    plot(ad,pXX,'+--','Color','#1f78b4',LineWidth=1.5,MarkerSize=5,HandleVisibility='off');
    xlabel('A-D $p$-value',Interpreter='latex',FontSize=15);
end
if logAxis == true
    set(gca, 'XScale', 'log');
    xline(0.005,'--',HandleVisibility='off');
    xline(0.1,'--',HandleVisibility='off');
end
ylim([0 200]);
set(gca,'YDir','reverse');
yticklabels([]);
grid minor;
legend(FontSize=15);

clr = 1:1:length(pXX);
subplot(1,5,[4.1 5])
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
hold on
plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
for i = 1:n2
    if strcmp(hypTest,'ks')
        if ks(2,i) < alphaKs
            plot(sk(i),ku(i),Marker="o",Color='k',HandleVisibility='off',MarkerSize=6);
        else
            plot(sk(i),ku(i),Marker="o",Color=[0.8 0.8 0.8],HandleVisibility='off');
        end
    else
        if ad(i) < alphaKs
            plot(sk(i),ku(i),Marker="o",Color='k',HandleVisibility='off',MarkerSize=6);
        else
            plot(sk(i),ku(i),Marker="o",Color=[0.8 0.8 0.8],HandleVisibility='off');
        end
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