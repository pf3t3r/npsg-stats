close all; clc; clear;

%% Load MaxMld and Chl-a (EpN)

epN = load("output\CTD\chla.mat").meanEpN;
epN = epN(1:61,:);
pIn = 0:2:120;
maxMld = load("mldVals.mat").maxMld;

% epNML = nan(61,329);
% 
% for i = 1:length(maxMld)
%     tmp = epN(1:(maxMld(i)/2),i);
%     epNML(1:(maxMld(i)/2),i) = tmp;
% end

% epNML = epN(1:(maxMld(1)/2),1);

%%
[ax,p,ks,obs,Sk,Ku,sd,rV,pV] = L1_helper_FLUORO(epN,pIn,maxMld);
