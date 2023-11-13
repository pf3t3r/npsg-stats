close all; clc; clear;
addpath("baroneRoutines\");

%% Load MaxMld and Chl-a (EpN)

% F2 = 131:329
epN = load("output\CTD\chla.mat").meanEpN;
epN = epN(1:76,131:329);
pIn = 0:2:150;
maxMld = load("mldVals.mat").maxMld;

%%
[ax,p,ks,obs,Sk,Ku,sd,rV,pV] = L1_helper_FLUORO(epN,pIn,maxMld);
sgtitle('CTD Chl a 01-21: L1');
exportgraphics(ax,'figures/L1/ctd/chla.png'); clear ax;
save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");

%% Load Temperature
ctdData = load("datafiles\ctd_iso_ALL.mat").ctd;
T = nan(329,76,31);
meanT = nan(76,329);

msng = [21, 48, 207, 218, 276];
tmp1 = 1:1:329;
tmp = setdiff(tmp1,msng);

for i = tmp
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

