%% Preamble
clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
set(0,'defaultAxesFontSize',10);

%% Need dates from CTD to match bottles to DCM
ctd = load('datafiles\ctd_iso_master2.mat').ctd;
datectd_HT = datetime([ctd.date],'ConvertFrom','datenum')';
datectd_UT = datectd_HT + hours(10);

%% Load and Clean Bottle Data

bottlePressure = importdata('data/hotbot-88_21.txt').data(:,4);
bottleChl = importdata('data/hotbot-88_21.txt').data(:,5);
botid = importdata('data/hotbot-88_21.txt').data(:,1);
timeA = importdata('data/hotbot-88_21.txt').data(:,2); % mmddyy

% Remove bottles taken at pressures below 2.5 db (within ~2.5 m of surface)
idRm = bottlePressure > 2.5;
bottlePressure = bottlePressure(idRm);
bottleChl = bottleChl(idRm);
botid = botid(idRm);
timeA = timeA(idRm);

% Remove bottles where chl-a concentration = 0
idZero = bottleChl == 0;
bottlePressure = bottlePressure(~idZero);
bottleChl = bottleChl(~idZero);
botid = botid(~idZero);
timeA = timeA(~idZero);

% Save cruise number (CRN) of each bottle - needed for Lagrangian analysis
tmp = num2str(botid);
bottleCRN = str2num(tmp(:,1:3));
clear tmp;

% Remove bottles from cruises 330 on (b/c fluorescence analysis not done)
for i = 1:9929
    if bottleCRN(i) > 329
        id329 = i - 1;
        break;
    end
end

bottlePressure = bottlePressure(1:id329);
bottleChl = bottleChl(1:id329);
botid = botid(1:id329);
timeA = timeA(1:id329);

clear idRm idZero id329 i;

%% time transformation

timeA(timeA==-9) = 0;
tTest2 = num2str(timeA,'%06d');
tTest3 = tTest2(:,5:6); %year
tTest3a = str2num(tTest3);
tTest4 = tTest2(:,3:4); %dd
tTest5 = tTest2(:,1:2); %mm
tTest6 = [];
tTest6 = [tTest3 tTest5 tTest4];
tTest7 = str2num(tTest6);
tTest8 = datetime(tTest7,'ConvertFrom','yyyymmdd');
tTest9 = num2str(19*ones(3256,1));
tTest10 = num2str(20*ones(9792-3256,1));
tTest11 = [tTest9; tTest10];
tTest12 = [tTest11 tTest3];
tTest13 = [tTest12 tTest5 tTest4];
tTest14 = str2num(tTest13);
tTest14(tTest14==19000000) = NaN;
tTest14(tTest14==20000000) = NaN;
tTest15 = datetime(tTest14,'ConvertFrom','yyyymmdd');
tIme = tTest15;

%% find cruise no.

ttBotIdString = num2str(botid);
ttCRN = str2num(ttBotIdString(:,1:3));
ttBottleCast = str2num(ttBotIdString(:,6:8)); % 100 => nan

for i = 1:329
    ttData2(i).cast = nan;
    ttData2(i).pL = nan;
end

for i = 1:329
    tmp = ttBottleCast(ttCRN==i);
    tmp2 = bottlePressure(ttCRN==i);
    tmp3 = bottleChl(ttCRN==i);    
    tmp2(tmp==100) = [];
    tmp3(tmp==100) = [];
    tmp(tmp==100) = []; % this order is important!
    tmp = unique(tmp);
    ttData2(i).cast = tmp;
    ttData2(i).pE = tmp2;
    ttData2(i).chla = tmp3;
    if isempty(ttData2(i).cast)
        ttData2(i).cast = nan;
    end
    if isempty(ttData2(i).chla)
        ttData2(i).chla = nan;
    end
    if isempty(ttData2(i).pE)
        ttData2(i).pE = nan;
    end
end

% Rewrite manually because here first cast = 3.
ttData2(34).cast = [9;10;14];

% for i = 1:329
%     tmp = ctd(i).pcm;
%     %disp(i);
%     if ~isnan(ttData2(i).cast)
%         ttData2(i).pcm = tmp(ttData2(i).cast);
%     end
%     if isempty(ttData2(i).pcm)
%         ttData2(i).pcm = nan;
%     end
%     if ~isnan(ttData2(i).pcm)
%         ttData2(i).pL = ttData2(i).pE - ttData2(i).pcm;
%         ttData2(i).pLbin = round(ttData2(i).pL,-1);
%     end
%     if isempty(ttData2(i).pLbin)
%         ttData2(i).pLbin = nan;
%     end
% end

tttest = [];
tttest2 = [];

% for i = 1:329
%     tmp = ttData2(i).pLbin;
%     tmp2 = ttData2(i).chla;
%     disp(i);
% %     for j = length() % continue later, fix pL
%     tttest = [tttest; tmp];
%     tttest2 = [tttest2; tmp2];
% end
% clear tmp tmp2;

% tttRange = unique(tttest);
% tttRange = tttRange(1:21);
% 
% for i = 1:length(tttRange)
%     tmp = tttest2(tttest==tttRange(i));
%     tmp(isnan(tmp)) = [];
%     [~,tttksVals(:,i),~] = statsplot2(tmp,'noplot');
% end

% axA = figure;
% 
% subplot(1,2,1)
% histogram(tttest);
% xlim([-150 150]);
% yline(100);
% set(gca,'XDir','reverse');
% set(gca,'view',[-90 90]);
% set(gca,'XAxisLocation','top');
% xlabel('Lagrangian Depth [dbar]');
% title('No. of Observations')

% subplot(1,2,2)
% plot(tttksVals(1,:),tttRange,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(tttksVals(2,:),tttRange,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(tttksVals(3,:),tttRange,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(tttksVals(4,:),tttRange,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% grid minor;
% ylim([-150 150]);
% set(gca,'YDir','reverse');
% legend('Location','best');
% xlabel('p-value');
% ylabel('Pressure [db]');
% title('KS Test');
% 
% sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle Chl-a');
% exportgraphics(axA,'figures/ks_NEWLAGBOT.png');
% clear axA;


% for i = 1:329
% end
% pLBin10 = round(pL,-1);

%% Bin Pressure

pBin5 = discretize(bottlePressure,2.5:5:202.5);
pBin10 = discretize(bottlePressure,2.5:10:202.5);
n = max(pBin5);
n10 = max(pBin10);

%% Visualise

% ax1 = figure;
% subplot(1,2,1)
% plot(bottleChl,bottlePressure,'r.','MarkerSize',4);
% xlabel('chl-a [ug/L]'); ylabel('Pressure [db]');
% set(gca,'YDir','reverse'); title('Raw Data');
% subplot(1,2,2)
% plot(bottleChl,pBin5,'b.','MarkerSize',4);
% xlabel('chl-a [ug/L]'); ylabel('Pressure [bin #]');
% set(gca,'YDir','reverse'); title('Binned Data');
% sgtitle('Fluorometric chl-a concentration (bottle data, 1988-2021)');
% exportgraphics(ax1,'figures/bottleChlaBinning.png');

%% Apply KS Test to chl-a across all pressures

ksEb5 = nan(5,40); ksEb10 = nan(5,20);
obsPerLevel_E5 = nan(40,1); obsPerLevel_E10 = nan(20,1);

for i = 1:n
    % find chl-a concentrations chla_i at binned pressure i
    chla_i = bottleChl(pBin5==i);
    % apply KS test to chla_i
    if length(chla_i) > 1
        [~,ksEb5(:,i),~] = statsplot2(chla_i,'noplot');
    end
    obsPerLevel_E5(i) = length(chla_i);
    clear chla_i;
end
clear n;

for i = 1:n10
    % find chl-a concentrations chla_i at binned pressure i
    chla_i10 = bottleChl(pBin10==i);
    % apply KS test to chla_i
    if length(chla_i10) > 1
        [~,ksEb10(:,i),~] = statsplot2(chla_i10,'noplot');
    end
    obsPerLevel_E10(i) = length(chla_i10);
    clear chla_i10;
end
clear n10;

depth5 = 5:5:200;
depth10 = 5:10:200;

ax2 = figure;
subplot(1,2,1)
barh(obsPerLevel_E5);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
set(gca,'YAxisLocation','right');
ylim([0 40]);
set(gca,"YTick",1:1:40,"YTickLabel",depth5)
title('No. of Observations');

for i = 1:40
    if obsPerLevel_E5(i) < 100
        ksEb5(:,i) = nan; % assign 0 instead of NaN b/c plot looks nicer :P
    end
end
for i = 1:20
    if obsPerLevel_E10(i) < 100
        ksEb10(:,i) = nan; % assign 0 instead of NaN b/c plot looks nicer :P
    end
end

% remove nan KS values for Eulerian 5dbar
ksEb5_f = ksEb5;
ksEb5_f = ksEb5_f(:,~all(isnan(ksEb5_f)));
depth5_f = depth5([1 2 3 5 9 12 15 17 19 20 21 23 25 27 30 35]);

ksEb5 = ksEb5_f;
% depth5 = depth5_f;
subplot(1,2,2)
plot(ksEb5(1,:),depth5_f,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEb5(2,:),depth5_f,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEb5(3,:),depth5_f,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEb5(4,:),depth5_f,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([0 200]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Eulerian Bottle chl-a [5 dbar bins]');
exportgraphics(ax2,'figures/ks_bottleEulerian5db.png');
clear ax2;

%%

% remove nan KS values for Eulerian 10dbar
ksEb10_f = ksEb10;
ksEb10_f = ksEb10_f(:,~all(isnan(ksEb10_f)));
depth10_f = depth10([1 2 3 4 5 6 8 9 10 11 12 13 14 15 18]);
ksEb10 = ksEb10_f;
% depth10 = depth10_f;

ax3 = figure;
subplot(1,2,1)
barh(obsPerLevel_E10);
hold on
xline(100);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
set(gca,'YAxisLocation','right');
ylim([0.5 20]);
set(gca,"YTick",1:1:20,"YTickLabel",depth10)
title('No. of Observations');

subplot(1,2,2)
plot(ksEb10(1,:),depth10_f,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEb10(2,:),depth10_f,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEb10(3,:),depth10_f,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEb10(4,:),depth10_f,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 195]);
grid minor;
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Eulerian Bottle chl-a [10 dbar bins]');
exportgraphics(ax3,'figures/ks_bottleEulerian10db.png');
clear ax3;

%% Lagrangian Bottle

% apply adjustment according to cruise avg lagrangian adjustments
% need to import this adjustment 
% for this I need to save the isoL and ctdL structs as mat files in
% nightTimeKs.m

ctdL = load("datafiles\lagrangianData.mat").ctdL;

% here I save the DCM for each cruise, but I need the DCM per cast...
for i = 1:329
    tmp = ctdL(i).dcmID;
    tmp(tmp==1) = nan;
    ctdL(i).dcmP = 2*tmp;
    dcmPcrn(i) = 2*round(mean(tmp,'omitnan')); %cruise-avg dcm
end

% bottlePressure - dcmP = lagP
% max no. of bottles per crn ~30 (?)
tCRN = 1; tBotId = 1; % initialise with CRN1, with starting ID = 1
for i = 2:9792
    if bottleCRN(i) > bottleCRN(i-1)
        tCRN = [tCRN string(bottleCRN(i))]; % cruise no. 
        tBotId = [tBotId i]; % starting ID of bottles that cruise
    end
end


dcmPcrn = dcmPcrn(str2double(tCRN));

tChlaByCrn = nan(112,319);
tPByCrn = nan(112,319);
tTime = NaT(112,319);
for i = 1:318
    tChlaByCrn(1:(tBotId(i+1)-tBotId(i)),i) = bottleChl(tBotId(i):(tBotId(i+1)-1));
    tPByCrn(1:(tBotId(i+1)-tBotId(i)),i) = bottlePressure(tBotId(i):(tBotId(i+1)-1));
    tTime(1:(tBotId(i+1)-tBotId(i)),i) = tIme(tBotId(i):(tBotId(i+1)-1));
end
tChlaByCrn(1:(9793-tBotId(319)),319) = bottleChl(tBotId(319):9792);
tPByCrn(1:(9793-tBotId(319)),319) = bottlePressure(tBotId(319):9792);
tTime(1:(9793-tBotId(319)),319) = tIme(tBotId(319):9792);

for i = 1:319
    pL(:,i) = tPByCrn(:,i) - dcmPcrn(i);
end

pLBin10 = round(pL,-1);

pL_Range = min(min(pLBin10)):10:max(max(pLBin10));

clear bottleCRN ctdL tBotId tCRN;

%% Hov

nb = 100;
ax4 = figure;
contourf(datenum(tTime),tPByCrn,tChlaByCrn,linspace(0,max(max(tChlaByCrn)),nb),'LineColor','auto');
set(gca,'YDir','reverse');
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Chlorophyll a [{\mu}g/L]';
xlabel('Time');
ylabel('Pressure [dbar]');
title('Bottle Chlorophyll a: Eulerian (1988-2021)');

exportgraphics(ax4,'figures/hov_botChla.png');

%%

pL_levels = nan(627,31); % pressure
chla_levels = nan(627,31); % chla
obsPerLevel_L = [];

for i = 1:31
    tmp = []; tmp2 = [];
    % find where pLBin10 == -160
    for j = 1:112
        for k = 1:319
            if pLBin10(j,k) == pL_Range(i)
                tmp = [tmp pLBin10(j,k)];
                tmp2 = [tmp2 tChlaByCrn(j,k)];
            end
        end
    end
    pL_levels(1:length(tmp),i) = tmp';
    chla_levels(1:length(tmp),i) = tmp2';
    obsPerLevel_L = [obsPerLevel_L length(tmp)];
    clear tmp tmp2;
end
clear j k;

%% KS it all

ksLb10 = getKS(chla_levels,31);

for i = 1:31
    if obsPerLevel_L(i) < 100
        ksLb10(:,i) = NaN; %again 0, not NaN, for nicer plot
    end
end
clear i;

%% figure this KS shit

pL_Range2 = pL_Range + 5;

ax5 = figure;

subplot(1,2,1)
histogram(pLBin10);
xlim([-150 150]);
yline(100);
set(gca,'XDir','reverse');
set(gca,'view',[-90 90]);
set(gca,'XAxisLocation','top');
xlabel('Lagrangian Depth [dbar]');
title('No. of Observations')

subplot(1,2,2)
plot(ksLb10(1,:),pL_Range2,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLb10(2,:),pL_Range2,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLb10(3,:),pL_Range2,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLb10(4,:),pL_Range2,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim([-150 150]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
ylabel('Pressure [db]');
title('KS Test');

sgtitle('Kolmogorov Smirnov Test: Lagrangian Bottle Chl-a');
exportgraphics(ax5,'figures/ksBotLag10db.png');
clear ax4 ax5;

%% Inspect for bimodality

% ax4 = figure;
% histfit(chla_i,[],'normal');
% hold on
% histfit(chla_i,[],'lognormal');
% histfit(chla_i,[],'weibull');
% histfit(chla_i,[],'gamma');
% legend();
% 
% exportgraphics(ax4,'figures/bottleLevelBimodality.png');