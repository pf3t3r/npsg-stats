% Most of this code is based on a faulty assumption, i.e. that the
% Kolomogorov-Smirnov Test may be applied to fluorescence given in
% isopycnal coordinates. But the isopycnals are defined only per cruise:
% between cruises, they change. Therefore the KS results are incorrect. So
% I will go thru this and remove this, and replace it with a proper
% approach whereby I find cruise average fluorescence on isopycnal
% coordinates but apply the KS test in pressure coordinates.

close all; clc; clear;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 29 15]);
set(0,'defaultAxesFontSize',11);
addpath("baroneRoutines\");

%%
% this is OK, just loading variables...

iso = load('datafiles\ctd_iso_master2.mat').iso;
ctd = load('datafiles\ctd_iso_master2.mat').ctd;

cruises = 329;
missingCruises = [21,207];
cruisesRecorded = linspace(1,cruises,cruises);
cruisesRecorded(missingCruises) = [];
clear cruises missingCruises;

%%
% 
% sig = [];
% f_iso = [];
% t_iso = [];
% pcm = [];
% for i = cruisesRecorded
%     tmpF = iso(i).f;
%     tmpT = ctd(i).date;
%     tmpPCM = ctd(i).pcm;
%     tmpS = iso(i).sig;
%     f_iso = [f_iso tmpF];
%     t_iso = [t_iso tmpT];
%     pcm = [pcm tmpPCM];
%     sig = [sig tmpS];
% end
% clear tmpF tmpT tmpPCM tmpS i;
% t_iso = datetime(t_iso,'ConvertFrom','datenum');
% 
% f_iso(f_iso < 0) = NaN;

%% test this

meanIsoF = [];
meanEulF = [];
for i = cruisesRecorded
    tmpI = mean([iso(i).f],2,"omitnan");
    meanIsoF = [meanIsoF tmpI];
    tmpE = mean([ctd(i).f],2,"omitnan");
    meanEulF = [meanEulF tmpE];
end
clear tmpI tmpE;

% Remove Negative Values
meanIsoF(meanIsoF<0) = nan;
meanEulF(meanEulF<0) = nan;

% This figure does not have casts removed: I need to incorporate the
% following sections
figure; plot(meanIsoF,[iso.sig]); set(gca,'YDir','reverse');
xlabel('Mean Fluorescence'); ylabel('Isopycnal Layer');
%% Lagrangian Conversion 
% this is incorrect I should not shift it this way

% first shift each individual cast
% inputs = offset, EulerianData, 

for i = cruisesRecorded
    % isopycnal
    tmp = iso(i).f(1:129,:);
    isoL(i).f = tmp;
    % pressure
    tmp2 = ctd(i).f(1:129,:);
    ctdL(i).f = tmp2;
end

for j = cruisesRecorded
    for i = 1:length(isoL(j).f(1,:)) % same length ctd and iso
        % iso
        tmp = isoL(j).f(:,i);
        [~,isoL(j).dcmID(i)] = max([isoL(j).f(:,i)]);
        isoL(j).offset(i) = 65 - isoL(j).dcmID(i);
        if isnan(isoL(j).offset(i))
            disp('i');
        else
            isoL(j).fL(:,i) = convertLagrangian(tmp,isoL(j).offset(i),329);
        end
        isoL(j).fLMean = mean(isoL(j).fL,2,'omitnan');

        % ctd
        tmp2 = ctdL(j).f(:,i);
        [~,ctdL(j).dcmID(i)] = max([ctdL(j).f(:,i)]);
        ctdL(j).offset(i) = 65 - ctdL(j).dcmID(i);
        if isnan(ctdL(j).offset(i))
            disp('i');
        else
            ctdL(j).fL(:,i) = convertLagrangian(tmp,ctdL(j).offset(i),329);
        end
        ctdL(j).fLMean = mean(ctdL(j).fL,2,'omitnan');
    end
end

%%

for i = cruisesRecorded
    meanFcm(i) = mean([ctd(i).fcm],2,'omitnan');
    meanPcm(i) = round(mean([ctd(i).pcm],2,'omitnan'));
end

% Step (1)
f_copy = meanEulF;
p_copy = meanPcm(cruisesRecorded);
colsToRemove = find(isnan(p_copy));

p_copy(colsToRemove) = [];
f_copy(:,colsToRemove) = [];

offset = 65 - round(p_copy/2);

tmp = f_copy(1:129,:);
meanLagF = nan(size(tmp));

% for i = 1:length(offset)
%     disp(i);
%     meanLagF(:,i) = circshift(tmp(:,i),offset(i));
%     if offset(i) > -1 && offset(i) < 40
%         meanLagF(1:offset(i),i) = NaN;
%     elseif offset(i) == -1
%         meanLagF(end,i) = NaN;
%     elseif offset(i) < -1 && offset(i) > -40
%         meanLagF((end+offset(i)):end,i) = NaN;
%     elseif abs(offset(i)) > 40
%         meanLagF(:,i) = NaN;
%     end
% end

%% Calculate KS
% Isopycnal Eulerian
for i = 1:129
    disp(i);
    tmp = meanIsoF(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksIsoF(:,i),~] = statsplot2(tmp,'noplot');
end

% Eulerian
for i = 1:129
    disp(i);
    tmp = meanEulF(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksEulF(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;


% Isopycnal Lagrangian
fLag = [isoL.fLMean];
fLag(fLag<0) = nan;
for i = 1:129
    disp(i);
    tmp = fLag(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksLagF(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
fLagP = [ctdL.fLMean];
fLagP(fLagP<0) = nan;
for i = 1:129
    disp(i);
    tmp = fLagP(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksLagFP(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;


%%

% The following figures show the KS test results for single cruise-averaged
% fluorescence profiles across time
% I have NOT removed the individual suspicious casts from the averaging, so
% 

ax = figure;

% Isopycnal Eulerian
subplot(1,4,1)
plot(ksIsoF(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksIsoF(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksIsoF(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksIsoF(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend('Location','best');
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Isopycnal Eulerian');

% Eulerian
subplot(1,4,2)
plot(ksEulF(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEulF(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEulF(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEulF(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian');

% Isopycnal Lagrangian
subplot(1,4,3)
plot(ksLagF(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLagF(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLagF(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLagF(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Isopycnal Lagrangian');

% Lagrangian
subplot(1,4,4)
plot(ksLagFP(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLagFP(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLagFP(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLagFP(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Lagrangian');

sgtitle('Kolmogorov Smirnov test: Cruise Average, 89-21');
exportgraphics(ax,'figures/ks_newInclIsopyc.png');

%% Night Time

lat = 22.75; lon = -158; 
rs = nan(329,31,2);
rsMean = nan(329,2);

for i = 1:329
    datetimeStruct(i).date = datetime(ctd(i).date,'ConvertFrom','datenum');
end

for i = 1:329
    for j = 1:length(datetimeStruct(i).date)
        [tmp,~,~,~,~,~] = suncycle(lat,lon,datetimeStruct(i).date(j));
        tmp = tmp - 10;
        for k = 1:2
            if tmp(k) < 0
                tmp(k) = tmp(k) + 24;
            end
        end
        rs(i,j,:) = tmp;
    end
    rsMean(i,1) = mean(rs(i,:,1),'omitnan');
    rsMean(i,2) = mean(rs(i,:,2),'omitnan');
end

rs2 = hours(rsMean);
rs2.Format = 'hh:mm';
% rs2 = rs2(cruisesRecorded,:);
% rs2(colsToRemove,:) = [];

%% cast at night or not

% reconvert to decimal hour for comparison
sunrise = hours(rs2(:,1));
sunset = hours(rs2(:,2));

for i = 1:329
    dayFrac = rem(datenum([ctd(i).date]),1);
    datetimeStruct(i).castTime = dayFrac*24;
end

datetimeStruct(i).nightID = [];
for i = cruisesRecorded
    for j = 1:length(datetimeStruct(i).date)
        if datetimeStruct(i).castTime(j) < sunrise(i)
            tmp = j;
            datetimeStruct(i).nightID = [datetimeStruct(i).nightID tmp];
        end
    end
end

meanIsoFN = [];
meanEulFN = [];
meanIsoLFN = [];
meanLFN = [];

for i = cruisesRecorded
    tmpIN = iso(i).f(:,datetimeStruct(i).nightID);
    tmpIN = mean(tmpIN,2,'omitnan');
    meanIsoFN = [meanIsoFN tmpIN];
    tmpEN = ctd(i).f(:,datetimeStruct(i).nightID);
    tmpEN = mean(tmpEN,2,"omitnan");
    meanEulFN = [meanEulFN tmpEN];
    tmpILN = isoL(i).fL(:,datetimeStruct(i).nightID);
    tmpILN = mean(tmpILN,2,'omitnan');
    meanIsoLFN = [meanIsoLFN tmpILN];
    tmpLN = ctdL(i).fL(:,datetimeStruct(i).nightID);
    tmpLN = mean(tmpLN,2,'omitnan');
    meanLFN = [meanLFN tmpLN];
end

meanIsoFN(meanIsoFN<0) = nan;
meanEulFN(meanEulFN<0) = nan;
meanIsoLFN(meanIsoLFN<0) = nan;
meanLFN(meanLFN<0) = nan;

%% NIGHT

% Isopycnal Eulerian
for i = 1:129
    disp(i);
    tmp = meanIsoFN(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksIsoFN(:,i),~] = statsplot2(tmp,'noplot');
end

% Eulerian
for i = 1:129
    disp(i);
    tmp = meanEulFN(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksEulFN(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;


% Isopycnal Lagrangian
for i = 1:129
    disp(i);
    tmp = meanIsoLFN(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksLagIsoFN(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:129
    disp(i);
    tmp = meanLFN(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksLagFN(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;


%% NIGHT
ax2 = figure;

% Isopycnal Eulerian
subplot(1,4,1)
plot(ksIsoFN(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksIsoFN(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksIsoFN(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksIsoFN(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend('Location','best');
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Isopycnal Eulerian');

% Eulerian
subplot(1,4,2)
plot(ksEulFN(1,:),[ctd(1).p(1:129)],'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksEulFN(2,:),[ctd(1).p(1:129)],'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksEulFN(3,:),[ctd(1).p(1:129)],'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksEulFN(4,:),[ctd(1).p(1:129)],'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Eulerian');

% Isopycnal Lagrangian
subplot(1,4,3)
plot(ksLagIsoFN(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLagIsoFN(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLagIsoFN(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLagIsoFN(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Isopycnal Lagrangian');

% Lagrangian
subplot(1,4,4)
plot(ksLagFN(1,:),[ctd(1).p(1:129)]-129,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksLagFN(2,:),[ctd(1).p(1:129)]-129,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksLagFN(3,:),[ctd(1).p(1:129)]-129,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksLagFN(4,:),[ctd(1).p(1:129)]-129,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca,'YDir','reverse');
xlabel('p-value');
ylim([-125 125]);
ylabel('Pressure [db]');
title('Lagrangian');

sgtitle('Kolmogorov Smirnov: Night Time Cruise Average, 89-21');
exportgraphics(ax2,'figures/ks_newInclIsopycNight.png');


%% MAYBE I should use the 'dodgy cast' removal part here ...
% 1. All those casts which have NaN as the pressure at the DCM, and
% 2. All those casts that were flagged for having anomalously large mean
% values (>0.5 ug/kg) and/or large single values (>2 ug/kg), and were
% confirmed visually as suspicious.

% Step (1)
% f_copy = f_iso;
% 
% p_copy = pcm;
% p_copy(isnan(pcm)) = [];
% 
% f_copy(:,isnan(pcm)) = [];
% f_isoR = f_copy(1:129,:);
% 
% t_copy = t_iso;
% t_copy(isnan(pcm)) = [];
% 
% % Step (2)
% dodgyCasts = [18 248 253 666 667 668 669 ...
%     1213 1214 1215 1216 1217 1218 1219 1220 1221 1222 1223 1224 ...
%     2124 2371 2434 2599 2610 2626 2898 2918 2926 2979 ...
%     3110 3133 3560 3561];
% f_isoR2 = f_isoR;
% f_isoR2(:,dodgyCasts) = [];
% t_copy(dodgyCasts) = [];

%% SECTIONS hereafter maybe don't use anymore

%% Find Lagrangian
% 
% [sig_cm,sig_cmID] = max(f_isoR2);
% f_iso_lag = NaN(129,length(f_isoR2(1,:)));
% offset_fISOL = 65 - sig_cmID;
% 
% for i = 1:length(f_isoR2(1,:))
%     f_iso_lag(:,i) = circshift(f_isoR2(1:129,i),offset_fISOL(i));
%     if offset_fISOL(i) > -1 && offset_fISOL(i) < 40
%         disp(i);
%         f_iso_lag(1:offset_fISOL(i),i) = NaN;
%     elseif offset_fISOL(i) == -1
%         f_iso_lag(end,i) = NaN;
%     elseif offset_fISOL(i) < -1 && offset_fISOL(i) > -40
%         disp(i);
%         f_iso_lag((end+offset_fISOL(i)):end,i) = NaN;
%     elseif abs(offset_fISOL(i)) > 40
%         f_iso_lag(:,i) = NaN;
%     end
% end

%% ALL CASTS
%% Eul + Lag All Cast: Calculate KS
% 
% % Eulerian
% for i = 1:129
%     disp(i);
%     tmp = f_isoR2(i,:);
%     tmp(isnan(tmp)) = [];
%     [~,ksISO(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;
% 
% % Lagrangian
% for i = 1:129
%     disp(i);
%     tmp = f_iso_lag(i,:);
%     tmp(isnan(tmp)) = [];
%     [~,ksISO_L(:,i),~] = statsplot2(tmp,'noplot');
% end
% clear tmp;

%% Eul + Lag All Cast: KS Figs

% incorrect
% sig_lag = sig(65) - sig(1:129);
% 
% ax = figure;
% ax.Position = [3 3 13 15];
% 
% % Eulerian
% subplot(1,2,1)
% plot(ksISO(1,:),sig(1:129),'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksISO(2,:),sig(1:129),'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksISO(3,:),sig(1:129),'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksISO(4,:),sig(1:129),'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% legend();
% xlim([0 0.8]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('\sigma^{\Theta} [kg m^{-3}]');
% title('Eulerian');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksISO_L(1,:),sig_lag,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksISO_L(2,:),sig_lag,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksISO_L(3,:),sig_lag,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksISO_L(4,:),sig_lag,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% % legend();
% xlim([0 0.6]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('\Delta \sigma^{\Theta} [kg m^{-3}]');
% title('Lagrangian');
% 
% sgtitle('KS on all casts (isopycnal, 89-21)');
% exportgraphics(ax,'figures/ks_iso_allCast.png');

%% NIGHT CASTS only
% %% Eul + Lag Night Cast: Calculate KS
% 
% nightCastID = load('datafiles\nightCast.mat').nightCastIDs;
% f_iso_en = f_isoR2(:,nightCastID);
% f_iso_ln = f_iso_lag(:,nightCastID);
% % t_n = t(nightCastIDs);
% 
% % Eulerian
% for i = 1:129
%     disp(i);
%     tmp = f_iso_en(i,:);
%     tmp(isnan(tmp)) = [];
%     if length(tmp) > 2
%         [~,ksIsoEN(:,i),~] = statsplot2(tmp,'noplot');
%     end
% end
% clear tmp;
% 
% % Lagrangian
% for i = 1:129
%     disp(i);
%     tmp = f_iso_ln(i,:);
%     tmp(isnan(tmp)) = [];
%     if length(tmp) > 2
%         [~,ksIsoLN(:,i),~] = statsplot2(tmp,'noplot');
%     end
% end
% clear tmp;

%% Eul + Lag Night Casts: KS Figs

% ax2 = figure;
% ax2.Position = [3 3 13 15];
% 
% % Eulerian
% subplot(1,2,1)
% plot(ksIsoEN(1,:),sig(1:129),'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksIsoEN(2,:),sig(1:129),'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksIsoEN(3,:),sig(1:129),'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksIsoEN(4,:),sig(1:129),'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% legend();
% xlim([0 0.8]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('\sigma^{\Theta} [kg m^{-3}]');
% title('Eulerian');
% 
% % Lagrangian
% subplot(1,2,2)
% plot(ksIsoLN(1,:),sig_lag,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ksIsoLN(2,:),sig_lag,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ksIsoLN(3,:),sig_lag,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ksIsoLN(4,:),sig_lag,'r.--','DisplayName','Gamma','MarkerSize',4);
% hold off
% % legend();
% xlim([0 0.6]);
% set(gca,'YDir','reverse');
% xlabel('p-value');
% ylabel('\Delta \sigma^{\Theta} [kg m^{-3}]');
% title('Lagrangian');
% 
% sgtitle('KS on night casts (isopycnal, 89-21)');
% exportgraphics(ax2,'figures/ks_iso_nightCast.png');

%% mean DCM
% 
% figure;
% plot(mean(f_isoR2,2,'omitnan'),sig(1:129),'Color',[0.6 0.6 0.6]);
% set(gca,'YDir','reverse');