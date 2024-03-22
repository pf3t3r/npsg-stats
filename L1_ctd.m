% Script to output L1 ctd results for the statistical analysis.

close all; clc; clear;
addpath("baroneRoutines\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 32 20]);

% Label for figures
% tmpT = "_NL";
% tmpT = "_S";
tmpT = "";

%% Load MaxMld and Chl-a (EpN)

% F2 = 131:329
epN = load("output\CTD\chla.mat").meanEpN;
epN = epN(1:76,131:329);
pIn = 0:2:150;
maxMld = load("mldVals.mat").maxMld;

%% load other variables

% Temperature T
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

figure
plot(meanT,pIn); set(gca,"YDir","reverse");
title("Temperature (C): 88-21");

figure
plot(meanSp,pIn); set(gca,"YDir","reverse");
title("Practical Salinity (g/kg): 88-21");

figure
plot(meanO2,pIn); set(gca,"YDir","reverse");
title("$O_2$ (mmol m$^{-3}$)",Interpreter="latex");

%% K-S

% CHL-A
ax = L1_helper_FLUORO(epN,pIn,maxMld);
sgtitle("CTD Chl a 01-21: L1");
exportgraphics(ax,"figures/L1/ctd/chla" + tmpT + ".png"); clear ax;
% save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");

% T
ax = L1_helper_FLUORO(meanT,pIn,maxMld);
sgtitle("CTD Temperature 88-21: L1");
exportgraphics(ax,"figures/L1/ctd/temp" + tmpT + ".png"); clear ax;
% save("output\L1\ctd\t.mat","p","ks","obs","Sk","Ku");

% SP
ax = L1_helper_FLUORO(meanSp,pIn,maxMld);
sgtitle("CTD S_P 88-21: L1");
exportgraphics(ax,"figures/L1/ctd/sp" + tmpT + ".png"); clear ax;
% save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");

% O2
ax = L1_helper_FLUORO(meanO2,pIn,maxMld);
sgtitle("CTD O2 88-21: L1");
exportgraphics(ax,"figures/L1/ctd/o2" + tmpT + ".png"); clear ax;
% save("output\L1\ctd\o2.mat","p","ks","obs","Sk","Ku");

%% A-D
tmpT = "-ad";

% CHL-A
ax = L1_helper_FLUORO(epN,pIn,maxMld,50,4,"ad");
sgtitle("CTD Chl a 01-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/ctd/chla" + tmpT + ".png"); clear ax;
% save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");

% T
ax = L1_helper_FLUORO(meanT,pIn,maxMld,50,4,"ad");
sgtitle("CTD Temperature 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/ctd/temp" + tmpT + ".png"); clear ax;
% save("output\L1\ctd\t.mat","p","ks","obs","Sk","Ku");

% SP
ax = L1_helper_FLUORO(meanSp,pIn,maxMld,50,4,"ad");
sgtitle("CTD S_P 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/ctd/sp" + tmpT + ".png"); clear ax;
% save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");

% O2
ax = L1_helper_FLUORO(meanO2,pIn,maxMld,50,4,"ad");
sgtitle("CTD O2 88-21: L1"+tmpT);
exportgraphics(ax,"figures/L1/ctd/o2" + tmpT + ".png"); clear ax;
% save("output\L1\ctd\o2.mat","p","ks","obs","Sk","Ku");

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
% [ax,p,ks,obs,Sk,Ku,rV,pV] = L1_helper_FLUORO(meanNO3,pIn,maxMld);
% sgtitle("CTD NO3- 88-21: L1");
% exportgraphics(ax,"figures/L1/ctd/notUsed/no3" + tmpT + ".png"); clear ax;
% % save("output\L1\ctd\no3.mat","p","ks","obs","Sk","Ku");