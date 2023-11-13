close all; clc; clear;
addpath("baroneRoutines\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 32 20]);

%% Load MaxMld and Chl-a (EpN)

% F2 = 131:329
epN = load("output\CTD\chla.mat").meanEpN;
epN = epN(1:76,131:329);
pIn = 0:2:150;
maxMld = load("mldVals.mat").maxMld;

%% Load Temperature
ctdData = load("datafiles\ctd_iso_ALL.mat").ctd;
msng = [21, 48, 207, 218, 276];
cR = 1:1:329;
cRm = setdiff(cR,msng);

T = nan(329,76,31);
meanT = nan(76,329);

for i = cRm
    tmp = ctdData(i).t(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            T(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    meanT(:,i) = mean(squeeze(T(i,:,:)),2,"omitnan");
end

%%
figure
plot(meanT,pIn); set(gca,'YDir','reverse');
title('Temperature (C): 88-21');

%% Load Salinity
% S missing data same as for T? Yes!

SP = nan(329,76,31);
meanSp = nan(76,329);

for i = cRm
    tmp = ctdData(i).sp(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            SP(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    meanSp(:,i) = mean(squeeze(SP(i,:,:)),2,"omitnan");
end

%%
figure
plot(meanSp,pIn); set(gca,'YDir','reverse');
title('Practical Salinity (g/kg): 88-21');

%% Load O2

O2 = nan(329,76,31);
meanO2 = nan(76,329);

for i = cRm
    tmp = ctdData(i).o(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            O2(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    meanO2(:,i) = mean(squeeze(O2(i,:,:)),2,"omitnan");
end

%% 
figure
plot(meanO2,pIn); set(gca,'YDir','reverse');
title('$O_2$ (mmol m$^{-3}$)',Interpreter='latex');

%% Load Nitrate (NO3-)

NO3 = nan(329,76,31);
meanNO3 = nan(76,329);

for i = cRm
    tmp = ctdData(i).n(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            NO3(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    meanNO3(:,i) = mean(squeeze(NO3(i,:,:)),2,"omitnan");
end

%%
figure
plot(meanNO3,pIn); set(gca,'YDir','reverse');
title('$NO_3^{-}$ (mmol m$^{-3}$)',Interpreter='latex');

%% L1 ANALYSIS

%% CHL-A
[ax,p,ks,obs,Sk,Ku,sd,rV,pV] = L1_helper_FLUORO(epN,pIn,maxMld);
sgtitle('CTD Chl a 01-21: L1');
% exportgraphics(ax,'figures/L1/ctd/chla.png'); clear ax;
save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");

%% T
[ax,p,ks,obs,Sk,Ku,sd,rV,pV] = L1_helper_FLUORO(meanT,pIn,maxMld);
sgtitle('CTD Temperature 88-21: L1');
% exportgraphics(ax,'figures/L1/ctd/temp.png'); clear ax;
save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");

%% O2
[ax,p,ks,obs,Sk,Ku,sd,rV,pV] = L1_helper_FLUORO(meanO2,pIn,maxMld);
sgtitle('CTD O2 88-21: L1');
% exportgraphics(ax,'figures/L1/ctd/o2.png'); clear ax;
save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");

%% NO3-
[ax,p,ks,obs,Sk,Ku,sd,rV,pV] = L1_helper_FLUORO(meanNO3,pIn,maxMld);
sgtitle('CTD NO3- 88-21: L1');
% exportgraphics(ax,'figures/L1/ctd/no3-.png'); clear ax;
save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");