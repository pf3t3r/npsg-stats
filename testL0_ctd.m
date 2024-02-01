clear; clc; close all;
addpath("baroneRoutines\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

%% CTD data

chlaCtd = load("output\CTD\chla.mat").meanEpN(1:101,131:329);
pCtd = 0:2:200;
tmpT = "";

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
        [~,ksC(:,i),~] = statsplot2(X_i,'noplot');
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
if max(kuC) > 10 & min(skC) < 0
    kurtLimB = max(kuC) + 1;
    skewLimA = min(skC) - 0.1;
    skewLimB = max(skC) + 0.1;
elseif max(kuC) > 10
    kurtLimB = max(kuC) + 1;
    skewLimB = max(skC) + 0.1;
elseif min(skC) < 0 
    skewLimA = min(skC) - 0.1;
elseif max(skC) > 2.5
    kurtLimB = max(kuC) + 1;
    skewLimB = max(skC) + 0.1;
else 
    kurtLimB = 10;
    skewLimA = 0;
    skewLimB = 2.5;
end

%%
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 18 15]);

% pXX = 5:10:195;
ax2 = figure;
subplot(1,2,1)
plot(ksC(2,:),pCtd,'+--','Color','#1f78b4',LineWidth=1.5,MarkerSize=5);
xline(0.05);
set(gca,'YDir','reverse');
ylabel('P [dbar]'); xlabel('K-S $p$-value','Interpreter','latex');

clr = 1:1:length(pCtd);
subplot(1,2,2)
plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
hold on
plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
scatter(skC,kuC,24,clr,"filled","o",HandleVisibility="off");
hold off
colormap(gca,cbrewer2("RdYlBu"));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:10:length(pCtd);
cbar.TickLabels = pCtd(1):20:pCtd(101);
cbar.Label.String = "P [dbar]";
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
ylabel('Kurtosis'); xlabel('Skewness');

sgtitle("L0: CTD Chl $a$","Interpreter","latex");
exportgraphics(ax2,"figures/L0/ctd/chla" + tmpT + ".png"); clear ax;

%%
% close all;

tmpAlpha = importdata('data/L0/alphaCaro_50-60.txt');
tmpBut19 = importdata('data\L0\but19_50-60.txt');
tmpHex19 = importdata('data\L0\hex19_50-60.txt');

alphaCar = tmpAlpha.data(:,5);
but19 = tmpBut19.data(:,5);
hex19 = tmpHex19.data(:,5);

%%
figure; histogram(alphaCar,max(alphaCar)); title('alpha-carotene');

figure; histogram(but19,max(but19)); title('But-19');

figure; histogram(hex19,max(hex19)); title('Hex-19');

%%
[mleA,ksA] = statsplot2(alphaCar);
[~,ksB] = statsplot2(but19);
[~,ksH] = statsplot2(hex19);


pd = makedist("Gamma","a",mleA(4,1),"b",mleA(4,2));

[~,pAdA] = adtest(alphaCar,"Distribution",pd);
[~,pAdB] = adtest(but19,"Distribution",pd);
[~,pAdH] = adtest(hex19,"Distribution",pd);

[rA,pA] = bbvuong(alphaCar);
[rB,pB] = bbvuong(but19);
[rH,pH] = bbvuong(hex19);

%%
dAcar = load("output\L1\acar.mat");
pOutA = dAcar.pOut;
cOutA = dAcar.cOut;

dHex19 = load("output\L1\hex19.mat");
pOutH = dHex19.pOut;
cOutH = dHex19.cOut;

dBut19 = load("output\L1\but19.mat");
pOutB = dBut19.pOut;
cOutB = dBut19.cOut;
%%
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 36 10]);

figure;
subplot(1,3,1)
histogram(alphaCar); title("alpha-caro");
subplot(1,3,2)
histogram(but19); title("but-19");
subplot(1,3,3)
histogram(hex19); title("hex-19");
sgtitle('50-60 dbar bin: L0');

%%
figure;
subplot(1,3,1)
histogram(cOutA(pOutA==6)); title("alpha-caro");
subplot(1,3,2)
histogram(cOutB(pOutB==6)); title("but-19");
subplot(1,3,3)
histogram(cOutH(pOutH==6)); title("hex-19");
sgtitle('50-60 dbar bin: L1');

%% remove outlier in hex-19 L0

hex19o = hex19;
hex19o(40) = [];
% show
figure
plot(hex19);
hold on
plot(hex19o);
hold off

[~,ksHo] = statsplot2(hex19o);
[rHo,pHo] = bbvuong(hex19o);