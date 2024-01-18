clear; clc; close all;
addpath("baroneRoutines\");
% addpath("data\L0\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

% Load 
tmp = importdata('data/L0/hplcChla_88-21_200.txt');
pIn = tmp.data(:,4);
chla = tmp.data(:,5);
idIn = tmp.data(:,1);

n = length(pIn);

tmpT = "";

chlaCtd = load("output\CTD\chla.mat").meanEpN(1:101,:);
pCtd = 0:2:200;
%% Bin
pB = discretize(pIn,0:10:200);
pX = nan(n,1);

for i = 1:n
    pX(i) = pB(i)*10 - 5;
end

%% K-S: bottle
threshold = 50;
n2 = 20;
ks = nan(5,n2);
sk = nan(1,n2); ku = nan(1,n2);
rV = nan(10,n2); pV = nan(10,n2);
obs = nan(1,n2);

% 4. Calculate KS p-value, skewness, kurtosis
for i = 1:n2
    % find concentration X_i at binned pressure i
    X_i = chla(pB==i);
    % apply KS test to X_i
    % change limit below to >3 to fix error with picoeu -> may change other
    % results
    if length(X_i) > 3
        [~,ks(:,i),~,~,~,~] = statsplot2(X_i,'noplot');
        [rV(:,i),pV(:,i)] = bbvuong(X_i);
        sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

%% K-S: CTD
threshold = 50;
n3 = 101;
ksC = nan(5,n3);
skC = nan(1,n3); kuC = nan(1,n3);
rVC = nan(10,n3); pVC = nan(10,n3);
obsC = nan(1,n3);

% 4. Calculate KS p-value, skewness, kurtosis
for i = 1:n3
    % find concentration X_i at binned pressure i
    X_i = chlaCtd(i,:);
    % apply KS test to X_i
    % change limit below to >3 to fix error with picoeu -> may change other
    % results
    X_i(isnan(X_i)) = [];
    if length(X_i) > 3
        [~,ksC(:,i),~,~,~,~] = statsplot2(X_i,'noplot');
        [rVC(:,i),pVC(:,i)] = bbvuong(X_i);
        skC(i) = skewness(X_i);
        kuC(i) = kurtosis(X_i);
    end
    obsC(i) = length(X_i);
    %clear X_i;
end


%%
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

%%

% ks(:,obs<threshold) = nan;

pXX = 5:10:195;
ax = figure;
subplot(1,3,1)
barh(obs);
xline(50); ylim([0.5 19.5]);
set(gca,"YDir","reverse");
title("No. of Observations");

subplot(1,3,2)
plot(ks(2,:),pXX,'+--','Color','#1f78b4',LineWidth=1.5,MarkerSize=5);
xline(0.05);
set(gca,'YDir','reverse');
ylabel('P [dbar]'); xlabel('p-value'); title('K-S');

clr = 1:1:length(pXX);
subplot(1,3,3)
plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
hold on
plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
scatter(sk,ku,24,clr,"filled","o",HandleVisibility="off");
hold off
colormap(gca,flipud(colormap("hot")));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:1:length(pXX);
cbar.TickLabels = pXX;
cbar.Label.String = "P [dbar]";
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
ylabel('Kurtosis'); xlabel('Skewness');
title("Skewness-Kurtosis");

sgtitle("L0: Bottle Chl-a");
exportgraphics(ax,"figures/L0/bottle/chla" + tmpT + ".png"); clear ax;

%%
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 18 15]);


% pXX = 5:10:195;
ax2 = figure;
subplot(1,2,1)
plot(ksC(2,:),pCtd,'+--','Color','#1f78b4',LineWidth=1.5,MarkerSize=5);
xline(0.05);
set(gca,'YDir','reverse');
ylabel('P [dbar]'); xlabel('$p$-value','Interpreter','latex'); title('K-S');

clr = 1:1:length(pCtd);
subplot(1,2,2)
plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
hold on
plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
scatter(skC,kuC,24,clr,"filled","o",HandleVisibility="off");
hold off
colormap(gca,flipud(colormap("hot")));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:10:length(pCtd);
cbar.TickLabels = pCtd(1):20:pCtd(101);
cbar.Label.String = "P [dbar]";
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
ylabel('Kurtosis'); xlabel('Skewness');
title("Skewness-Kurtosis");

sgtitle("L0: CTD Chl-a");
exportgraphics(ax2,"figures/L0/ctd/chla" + tmpT + ".png"); clear ax;