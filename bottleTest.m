%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 15 15]);
set(0,'defaultAxesFontSize',10);

%% Load bottle data
chlBotData = importdata('data/hotbot-88_21.txt').data;
bottlePressure = chlBotData(:,4);
bottleChl = chlBotData(:,5);

%% Bin fluorometric chl-a in 5 db intervals
% ... such that central value of each bin is 5,10,15,etc.
botbot = bottlePressure(bottlePressure>2.4);

binnedPressure = discretize(botbot,2.5:5:202.5);

% Find no. of measurements at each depth
histSet = zeros(40,1);
x = 5:5:200;
n = length(botbot);
for i = 1:40
    tmp = 0;
    for j = 1:n
        if binnedPressure(j) == i
            tmp = tmp+1;
        end
    end
    histSet(i) = tmp;
    clear tmp;
end

figure;
barh(x,histSet,'HandleVisibility','off');
xline(100,'r-','DisplayName','Threshold');
set(gca,'YDir','reverse');
legend('Location','best');
title('Fluorometric Chl-a: no. of observations at each depth class (Eulerian)');

%% Separate bottle measurements into distinct 'casts'

% % find start and end of bottle 'cast'
% botcaststart = 1;
% botcastend = [];
% 
% for i=1:n-1
%     if bottlePressure(i+1) - bottlePressure(i) < 0
%         botcastend = [botcastend i];
%         botcaststart = [botcaststart i+1];
%     end
% end
% 
% botcastend = [botcastend n];
% 
% castLength = [];
% for i=1:1662
%     castLength = [castLength 1+botcastend(i)-botcaststart(i)];
% end
% maxCast = max(castLength);

%% plot above
% 
% % Save the casts
% chlBotCasts = nan([maxCast,1662]);
% pBotCasts = nan([maxCast,1662]);
% chlBotCasts(1:length(bottleChl(1:4)),1) = bottleChl(1:4);
% pBotCasts(1:length(bottleChl(1:4)),1) = bottlePressure(1:4);
% 
% for i=2:1662
%     chlBotCasts(1:length(bottleChl(botcaststart(i):botcastend(i))),i) = bottleChl(botcaststart(i):botcastend(i));
%     pBotCasts(1:length(bottleChl(botcaststart(i):botcastend(i))),i) = bottlePressure(botcaststart(i):botcastend(i));
% end
% hold off
% 
% valbot = NaN(1662,1);
% idbot = NaN(1662,1);
% chlBotCastsID = [];
% pBotCastsID = [];
% 
% % find DCM
% for i=1:1662
%     if nnz(~isnan((chlBotCasts(:,i)))) > 3
%         [valbot(i),idbot(i)] = max(chlBotCasts(:,i));
%         chlBotCastsID = [chlBotCastsID i];
%         pBotCastsID = [pBotCastsID i];
%     end
% end

%% Select only 'casts' with 4+ measurements

% chlBotCasts2 = chlBotCasts(:,chlBotCastsID);
% pBotCasts2 = pBotCasts(:,pBotCastsID);
% 
% valbot2 = valbot(chlBotCastsID);
% idbot2 = idbot(chlBotCastsID);
%%
% figure;
% plot(bottleChl(1:4),bottlePressure(1:4),'Color',[0.6 0.6 0.6]);
% set(gca,'YDir','reverse');
% hold on
% for i = 1000:1100
%     plot(bottleChl(botcaststart(i):botcastend(i)),bottlePressure(botcaststart(i):botcastend(i)),'Color',[0.6 0.6 0.6],'LineWidth',0.1);
%     plot(chlBotCasts(idbot(i),i),pBotCasts(idbot(i),i),'.r');
% end
% hold off

% %% 
% figure;
% plot(chlBotCasts2(1:4,1),pBotCasts2(1:4,1),'Color',[0.6 0.6 0.6]);
% set(gca,'YDir','reverse');
% hold on
% for i = 1:length(chlBotCastsID)
%     plot(chlBotCasts2(:,i),pBotCasts2(:,i),'Color',[0.6 0.6 0.6]);
% %     plot(chlBotCasts2(idbot2(i),i),pBotCasts2(idbot2(i),i),'.r');
% end
% plot(mean(chlBotCasts2,2,"omitnan"),mean(pBotCasts2,2,"omitnan"),'.r');
% hold off

%% Bin the separated casts
% pBinned = discretize(pBotCasts2,2.5:5:202.5);

%% sort the pressure

[B,I] = sort(binnedPressure);
sortedChl = bottleChl(I);

%%
% figure;
% plot(sortedChl,B,'r.');
% set(gca,'YDir','reverse');

%% find start/end of level
endOfLevel = [0];
for i = 1:9939-1
    if B(i+1) > B(i)
        endOfLevel = [endOfLevel i];
    end
end

%% KS for Eulerian

n = length(histSet);

ksEb = zeros(5,n);
for i = 1:39
    disp(i);
    tmp = sortedChl(endOfLevel(i)+1:endOfLevel(i+1));
    tmp(isnan(tmp) | tmp<=0) = [];
    [~,ksEb(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;

%% KS figure

% Eulerian
% figure;
figure;
plot(ksEb(1,:),x,'Color',[0 0 0],'DisplayName','Normal');
hold on
plot(ksEb(2,:),x,'Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal');
plot(ksEb(3,:),x,'Color','red','DisplayName','Weibull');
plot(ksEb(4,:),x,'Color','red','LineStyle','--','DisplayName','Gamma');
% plot(ksEb(5,:),x,'Color','red','LineStyle',':','DisplayName','Exp');
hold off
set(gca,'YDir','reverse');
legend();
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian KS-Test');