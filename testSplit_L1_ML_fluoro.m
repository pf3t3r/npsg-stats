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

%% Load other parameters

%% Temperature
ctdData = load("datafiles\ctd_iso_ALL.mat").ctd;
T = nan(329,501,31);

for i = 1:1
    tmp = ctdData(i).t;
    for j = 1:length(tmp(1,:))
        T(1,:,j) = tmp(:,j);
    end
end

meanT = mean(squeeze(T(1,:,:)),2,"omitnan");