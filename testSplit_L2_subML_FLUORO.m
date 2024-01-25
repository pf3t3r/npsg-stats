clear; clc; close all; addpath("baroneRoutines\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

%% Load MLD and DCM
maxMld = load('mldVals.mat').maxMld;
dcm = load("dcm.mat").meanPcm;  

%% Assign lower pressure we measure until...

lowerP = 129;
pIn = 0:2:2*(lowerP-1);

%% Chlorophyll a: 88-21

% 1. Load data
chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,131:329);

% Extract data below ML and centre on DCM; calculate KS p-value, skewness,
% and kurtosis; and plot.
[ax,pL,ks,obs,sk,ku,pV,rV,tr2] = L2_helper_FLUORO(chla,pIn,maxMld,dcm);
sgtitle('[Chl a] 01-21: L2');
exportgraphics(ax,'figures/L2/ctd/chla.png');
save("output\L2\ctd\chla.mat","pL","ks","obs","sk","ku","pV","rV");

% ax2 = figure;
% plot(chla,pL,LineStyle=":",Color=[0.6 0.6 0.6]);
% set(gca,"YDir","reverse");
% xlabel('Chl-a [mg m^{-3}]'); ylabel('P [dbar]'); title('Chl-a','P = 0 dbar => DCM');
% exportgraphics(ax2,'figures/L2/chla.png');

clearvars -except maxMld dcm lowerP pIn chla sk ku tr2;

%% Temperature: 88-21

% Load Data
ctdData = load("datafiles\ctd_iso_ALL.mat").iso;
msng = [21, 48, 207, 218, 276];
cR = 1:1:329;
cRm = setdiff(cR,msng);
clear msng cR;

Tin = nan(329,lowerP,31);
T = nan(lowerP,329);

for i = cRm
    tmp = ctdData(i).t(1:lowerP,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            Tin(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    T(:,i) = mean(squeeze(Tin(i,:,:)),2,"omitnan");
end

[ax,pL,ks,obs,sk,ku,pV,rV] = L2_helper_FLUORO(T,pIn,maxMld,dcm);
sgtitle('T 88-21: L2');
exportgraphics(ax,'figures/L2/ctd/T.png'); clear ax;
save("output\L2\ctd\T.mat","pL","ks","obs","sk","ku","pV","rV");

% ax2 = figure;
% plot(T,pL,LineStyle=":",Color=[0.6 0.6 0.6]);
% set(gca,"YDir","reverse");
% xlabel('T [Celsius]'); ylabel('P [dbar]'); title('Temperature','P = 0 dbar => DCM');
% exportgraphics(ax2,'figures/L2/T.png');

clearvars -except maxMld dcm lowerP pIn chla T ctdData cRm;

%% Salinity: 88-21

SPin = nan(329,lowerP,31);
Sp = nan(lowerP,329);

for i = cRm
    tmp = ctdData(i).sp(1:lowerP,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            SPin(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    Sp(:,i) = mean(squeeze(SPin(i,:,:)),2,"omitnan");
end

[ax,pL,ks,obs,sk,ku,pV,rV] = L2_helper_FLUORO(Sp,pIn,maxMld,dcm);
sgtitle('S_p 88-21: L2');
exportgraphics(ax,'figures/L2/ctd/Sp.png'); clear ax;
save("output\L2\ctd\Sp.mat","pL","ks","obs","sk","ku","pV","rV");

% ax2 = figure;
% plot(Sp,pL,LineStyle=":",Color=[0.6 0.6 0.6]);
% set(gca,"YDir","reverse");
% xlabel('S_p [g/kg]'); ylabel('P [dbar]'); title('Salinity','P = 0 dbar => DCM');
% exportgraphics(ax2,'figures/L2/Sp.png');

% clearvars -except maxMld dcm lowerP pIn chla T Sp ctdData cRm;

%% O2: 88-21

O2in = nan(329,lowerP,31);
o2 = nan(lowerP,329);

for i = cRm
    tmp = ctdData(i).o(1:lowerP,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            O2in(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    o2(:,i) = mean(squeeze(O2in(i,:,:)),2,"omitnan");
end

[ax,pL,ks,obs,sk,ku,pV,rV] = L2_helper_FLUORO(o2,pIn,maxMld,dcm);
sgtitle('O_2 88-21: L2');
exportgraphics(ax,'figures/L2/ctd/o2.png');
save("output\L2\ctd\o2.mat","pL","ks","obs","sk","ku","pV","rV");

% ax2 = figure;
% plot(o2,pL,LineStyle=":",Color=[0.6 0.6 0.6]);
% set(gca,"YDir","reverse");
% xlabel('O_2 [mmol m^{-3}]'); ylabel('P [dbar]'); title('O_2','P = 0 dbar => DCM');
% exportgraphics(ax2,'figures/L2/o2.png');

clearvars -except maxMld dcm lowerP pIn chla T Sp o2 ctdData cRm;

%% No3-: 88-21

NO3in = nan(329,lowerP,31);
no3 = nan(lowerP,329);

for i = cRm
    tmp = ctdData(i).n(1:lowerP,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            NO3in(i,:,j) = tmp(:,j);
        end
    end
end

for i = 1:329
    no3(:,i) = mean(squeeze(NO3in(i,:,:)),2,"omitnan");
end

[ax,pL,ks,obs,sk,ku,pV,rV] = L2_helper_FLUORO(no3,pIn,maxMld,dcm);
sgtitle('NO_3^{-} 88-21: L2');
exportgraphics(ax,'figures/L2/ctd/no3.png');
save("output\L2\ctd\no3.mat","pL","ks","obs","sk","ku","pV","rV");

% ax2 = figure;
% plot(no3,pL,LineStyle=":",Color=[0.6 0.6 0.6]);
% set(gca,"YDir","reverse");
% xlabel('NO_3^{-} [mmol m^{-3}]'); ylabel('P [dbar]'); title('NO_3','P = 0 dbar => DCM');
% exportgraphics(ax2,'figures/L2/no3.png');

clearvars -except maxMld dcm lowerP pIn chla T Sp o2 no3 ctdData cRm;