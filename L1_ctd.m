% Script to output L1 ctd results for the statistical analysis.

close all; clc; clear;
addpath("baroneRoutines\"); addpath("func\"); addpath("output\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 32 20]);

% Test Cases
principleAnalysis = false;       % main analysis
seasonalAnalysis = true;       % seasonality of statistics
logAxes = true;                 % output p-values as log values (true)
if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Load MaxMld and Chl-a (EpN) and CTD data

% F2 = 131:329
epN = load("output\CTD\chla.mat").meanEpN;
epN = epN(1:76,131:329);
pIn = 0:2:150;
maxMld = load("mldVals.mat").maxMld;

ctdData = load("datafiles\ctd_iso_ALL.mat").ctd;
msng = [21, 48, 207, 218, 276];
cR = 1:1:329;
cRm = setdiff(cR,msng);

% Temperature T
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

% Salinity SP
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

% O2
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

%% Show vertical profiles

% figure
% plot(meanT,pIn); set(gca,"YDir","reverse");
% title("Temperature (C): 88-21");
% 
% figure
% plot(meanSp,pIn); set(gca,"YDir","reverse");
% title("Practical Salinity (g/kg): 88-21");
% 
% figure
% plot(meanO2,pIn); set(gca,"YDir","reverse");
% title("$O_2$ (mmol m$^{-3}$)",Interpreter="latex");

%% Principal Analysis
if principleAnalysis == true
    %  K-S
    tmpT = "";
    
    % CHL-A
    ax = L1_ctdHelper(epN,pIn,maxMld);
    sgtitle("CTD Chl a 01-21: L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");
    
    % T
    [ax,mldCon,rV,pV,ks,vuongRes] = L1_ctdHelper(meanT,pIn,maxMld);
    sgtitle("CTD Temperature 88-21: L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\t.mat","p","ks","obs","Sk","Ku");
    
    % SP
    ax = L1_ctdHelper(meanSp,pIn,maxMld);
    sgtitle("CTD S_P 88-21: L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"sp" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");
    
    % O2
    [ax,mldCon] = L1_ctdHelper(meanO2,pIn,maxMld);
    sgtitle("CTD O2 88-21: L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"o2" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\o2.mat","p","ks","obs","Sk","Ku");
    
    % A-D
    tmpT = "-ad";
    
    % CHL-A
    ax = L1_ctdHelper(epN,pIn,maxMld,50,4,"ad");
    sgtitle("CTD Chl a 01-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");
    
    % T
    ax = L1_ctdHelper(meanT,pIn,maxMld,50,4,"ad");
    sgtitle("CTD Temperature 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\t.mat","p","ks","obs","Sk","Ku");
    
    % SP
    [ax,mldCon,rV,pV,ad,V] = L1_ctdHelper(meanSp,pIn,maxMld,50,4,"ad");
    sgtitle("CTD S_P 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"sp" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");
    
    % O2
    ax = L1_ctdHelper(meanO2,pIn,maxMld,50,4,"ad");
    sgtitle("CTD O2 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"o2" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\o2.mat","p","ks","obs","Sk","Ku");
end

%% Set-up seasonal test.

% find "average" month of a cruise, then give an ID to that cruise saying
% which season it is (Spring, Summer, Autumn, or Winter).

avgMonth = nan(329,1);
winter = nan(329,1);
spring = nan(329,1);
summer = nan(329,1);
autumn = nan(329,1);

for i = 1:329
    tmp = round(mean(month(datetime(ctdData(i).date,"ConvertFrom","datenum"))));
    
    if (tmp == 12) || (tmp == 1) || (tmp == 2)
        winter(i) = 1;
    end
    if (tmp == 3) || (tmp == 4) || (tmp == 5)
        spring(i) = 1;
    end
    if (tmp == 6) || (tmp == 7) || (tmp == 8)
        summer(i) = 1;
    end
    if (tmp == 9) || (tmp == 10) || (tmp == 11)
        autumn(i) = 1;
    end

    avgMonth(i) = tmp;
end

% ASTRO seasons -> winter = jfm; spring = amj; summer = jas; autumn = ond.

winIds = []; sprIds = []; sumIds = []; autIds = []; 
for i = 1:329
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

% Calculate average mld per season
mldW = mean(maxMld(winIds),"omitnan");
mldSp = mean(maxMld(sprIds),"omitnan");
mldSu = mean(maxMld(sumIds),"omitnan");
mldA = mean(maxMld(autIds),"omitnan");

% Average dcm per season
dcm = load("dcm.mat").meanPcm;  
dcmW = mean(dcm(winIds),"omitnan");
dcmSp = mean(dcm(sprIds),"omitnan");
dcmSu = mean(dcm(sumIds),"omitnan");
dcmA = mean(dcm(autIds),"omitnan");

%% Test Seasonal Analysis

% Naming scheme
% -01 = winter, -02 = spring, -03 = summer, -04 = autumn
adThresh = 30; % only for chl-a in A-D

if seasonalAnalysis == true

    % K-S
    % winter
    tmpT = "-01";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,winIds(36:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,50,4,"ks");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,winIds),pIn,maxMld,50,4,"ks");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,winIds),pIn,maxMld,50,4,"ks");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,winIds),pIn,maxMld,50,4,"ks");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % spring
    tmpT = "-02";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sprIds(33:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,50,4,"ks");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,sprIds),pIn,maxMld,50,4,"ks");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,sprIds),pIn,maxMld,50,4,"ks");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,sprIds),pIn,maxMld,50,4,"ks");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % summer
    tmpT = "-03";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(33:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,50,4,"ks");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,sumIds),pIn,maxMld,50,4,"ks");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,sumIds),pIn,maxMld,50,4,"ks");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,sumIds),pIn,maxMld,50,4,"ks");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % autumn
    tmpT = "-04";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(30:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,50,4,"ks");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,autIds),pIn,maxMld,50,4,"ks");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,autIds),pIn,maxMld,50,4,"ks");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,autIds),pIn,maxMld,50,4,"ks");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % A-D
    % winter
    tmpT = "-ad-01";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,winIds(36:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,adThresh,4,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,winIds),pIn,maxMld,50,4,"ad");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,winIds),pIn,maxMld,50,4,"ad");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,winIds),pIn,maxMld,50,4,"ad");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % spring
    tmpT = "-ad-02";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sprIds(33:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,adThresh,4,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,sprIds),pIn,maxMld,50,4,"ad");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,sprIds),pIn,maxMld,50,4,"ad");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,sprIds),pIn,maxMld,50,4,"ad");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % summer
    tmpT = "-ad-03";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(33:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,adThresh,4,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,sumIds),pIn,maxMld,50,4,"ad");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,sumIds),pIn,maxMld,50,4,"ad");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,sumIds),pIn,maxMld,50,4,"ad");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % autumn
    tmpT = "-ad-04";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(30:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,adThresh,4,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,autIds),pIn,maxMld,50,4,"ad");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,autIds),pIn,maxMld,50,4,"ad");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,autIds),pIn,maxMld,50,4,"ad");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");
end


%% Unused.
% % Load Nitrate (NO3-)
% 
% NO3 = nan(329,76,31);
% meanNO3 = nan(76,329);
% 
% for i = cRm
%     tmp = ctdData(i).n(1:76,:);
%     if length(tmp) > 3
%         for j = 1:length(tmp(1,:))
%             NO3(i,:,j) = tmp(:,j);
%         end
%     end
% end
% 
% for i = 1:329
%     meanNO3(:,i) = mean(squeeze(NO3(i,:,:)),2,"omitnan");
% end
% %
% figure
% plot(meanNO3,pIn); set(gca,"YDir","reverse");
% title("$NO_3^{-}$ (mmol m$^{-3}$)",Interpreter="latex");
% %% NO3-
% [ax,p,ks,obs,Sk,Ku,rV,pV] = L1_ctdHelper(meanNO3,pIn,maxMld);
% sgtitle("CTD NO3- 88-21: L1");
% exportgraphics(ax,"figures/L1/ctd/notUsed/no3" + tmpT + ".png"); clear ax;
% % save("output\L1\ctd\no3.mat","p","ks","obs","Sk","Ku");