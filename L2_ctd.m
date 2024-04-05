% Script to output L2 ctd results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\"); addpath("output\");
set(groot,'defaultFigureUnits','centimeters','defaultFigurePosition',[3 3 28 15]);

% Load MLD and DCM
maxMld = load('mldVals.mat').maxMld;
dcm = load("dcm.mat").meanPcm;  

% Assign lower bound on pressure
lowerP = 129;
pIn = 0:2:2*(lowerP-1);

% Cases studied here
principalAnalysis = false;
seasonalAnalysis = true;
logAxes = true;                 % output p-values as log values (true)
if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Load hydrographical variables

ctdData = load("datafiles\ctd_iso_ALL.mat").iso;
timeData = load("datafiles\ctd_iso_ALL.mat").ctd;
msng = [21, 48, 207, 218, 276];
cR = 1:1:329;
cRm = setdiff(cR,msng);
clear msng cR;

% Temperature
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

% Salinity
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

clear i j cRm O2in SPin Tin tmp;
% keep dcm lowerP maxMld pIn o2 Sp T



%% Seasonal Analysis

if seasonalAnalysis == true

    % Set-up seasonal test.

    % find "average" month of a cruise, then give an ID to that cruise 
    % saying which season it is (Spring, Summer, Autumn, or Winter).
    
    avgMonth = nan(329,1);
    winter = nan(329,1);
    spring = nan(329,1);
    summer = nan(329,1);
    autumn = nan(329,1);
    
    for i = 1:329
        tmp = round(mean(month(datetime(timeData(i).date,"ConvertFrom","datenum"))));
        
        if (tmp == 1) || (tmp == 2) || (tmp == 3)
            winter(i) = 1;
        end
        if (tmp == 4) || (tmp == 5) || (tmp == 6)
            spring(i) = 1;
        end
        if (tmp == 7) || (tmp == 8) || (tmp == 9)
            summer(i) = 1;
        end
        if (tmp == 10) || (tmp == 11) || (tmp == 12)
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

    %% Apply Seasonal Analysis

    % K-S
    % winter
    tmpT = "-01";

    % chla
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,winIds(36:end));
    [ax,kLim,tr2,ks] = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ks",[-60 60],25);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L2/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L2_ctdHelper(T(:,winIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png"); clear ax;

    % Sp
    ax = L2_ctdHelper(Sp(:,winIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png"); clear ax;

    % O2
    ax = L2_ctdHelper(o2(:,winIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");

    % spring
    tmpT = "-02";

    % chla
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,sprIds(33:end));
    ax = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ks",[-60 60],25);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L2/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L2_ctdHelper(T(:,sprIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png");

    % Sp
    ax = L2_ctdHelper(Sp(:,sprIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png");

    % O2
    ax = L2_ctdHelper(o2(:,sprIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");

    % summer
    tmpT = "-03";

    % chla
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,sumIds(33:end));
    ax = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ks",[-60 60],25);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L2/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L2_ctdHelper(T(:,sumIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png");

    % Sp
    ax = L2_ctdHelper(Sp(:,sumIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png");

    % O2
    ax = L2_ctdHelper(o2(:,sumIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");

    % autumn
    tmpT = "-04";

    % chla
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,autIds(30:end));
    ax = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ks",[-60 60],25);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L2/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L2_ctdHelper(T(:,autIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png");

    % Sp
    ax = L2_ctdHelper(Sp(:,autIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png");

    % O2
    ax = L2_ctdHelper(o2(:,autIds),pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");

    %%%

    % A-D
    % winter
    tmpT = "-ad-01";

    % chla
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,winIds(36:end));
    ax = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ad",[-60 60],25);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L2/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L2_ctdHelper(T(:,winIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png"); clear ax;

    % Sp
    ax = L2_ctdHelper(Sp(:,winIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png"); clear ax;

    % O2
    ax = L2_ctdHelper(o2(:,winIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");

    % spring
    tmpT = "-ad-02";

    % chla
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,sprIds(33:end));
    ax = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ad",[-60 60],25);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L2/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L2_ctdHelper(T(:,sprIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png");

    % Sp
    ax = L2_ctdHelper(Sp(:,sprIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png");

    % O2
    ax = L2_ctdHelper(o2(:,sprIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");

    % summer
    tmpT = "-ad-03";

    % chla
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,sumIds(33:end));
    ax = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ad",[-60 60],25);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L2/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L2_ctdHelper(T(:,sumIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png");

    % Sp
    ax = L2_ctdHelper(Sp(:,sumIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png");

    % O2
    ax = L2_ctdHelper(o2(:,sumIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");

    % autumn
    tmpT = "-ad-04";

    % chla
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,autIds(30:end));
    ax = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ad",[-60 60],25);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L2/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L2_ctdHelper(T(:,autIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png");

    % Sp
    ax = L2_ctdHelper(Sp(:,autIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png");

    % O2
    ax = L2_ctdHelper(o2(:,autIds),pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");

end

%% Principal Analysis
if principalAnalysis == true

    % K-S
    tmpT = "";
    
    % Chlorophyll a (88-21)
    chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,131:329);
    %[ax,pL,ks,obs,sk,ku,pV,rV,~,~,tr2] = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ks",[-60 60]);
    [ax,kLim,tr2,ks] = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('[Chl a] 01-21: L2');
    exportgraphics(ax,"figures/L2/ctd/log/chla" + tmpT + ".png");
    save("output\L2\ctd\chla.mat","pL","ks","obs","sk","ku","pV","rV");
    clear ax ks ku obs pL pV rV sk;
    
    % Temperature (88-21)
    [ax,pL,ks,obs,sk,ku,pV,rV] = L2_ctdHelper(T,pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('T 88-21: L2');
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png"); clear ax;
    save("output\L2\ctd\T.mat","pL","ks","obs","sk","ku","pV","rV");
    clear ax ks ku obs pL pV rV sk;
    
    % Practical Salinity (88-21)
    [ax,pL,ks,obs,sk,ku,pV,rV] = L2_ctdHelper(Sp,pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('S_p 88-21: L2');
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png"); clear ax;
    save("output\L2\ctd\Sp.mat","pL","ks","obs","sk","ku","pV","rV");
    clear ax ks ku obs pL pV rV sk;
    
    % O2 (88-21)
    [ax,pL,ks,obs,sk,ku,pV,rV] = L2_ctdHelper(o2,pIn,maxMld,dcm,4,"ks",[-60 60]);
    sgtitle('O_2 88-21: L2');
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");
    save("output\L2\ctd\o2.mat","pL","ks","obs","sk","ku","pV","rV");
    clear ax ks ku obs pL pV rV sk;
    
    % A-D
    tmpT = "-ad";
    
    % Chlorophyll a (88-21)
    % chla = load("output\CTD\chla.mat").meanLiN(1:lowerP,131:329);
    [ax,pL,ks,obs,sk,ku,pV,rV] = L2_ctdHelper(chla,pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('[Chl a] 01-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/chla" + tmpT + ".png");
    save("output\L2\ctd\chla.mat","pL","ks","obs","sk","ku","pV","rV");
    clear ax lowerP ks ku obs pL pV rV sk;
    
    % Temperature (88-21)
    [ax,pL,ks,obs,sk,ku,pV,rV] = L2_ctdHelper(T,pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('T 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/T" + tmpT + ".png");
    save("output\L2\ctd\T.mat","pL","ks","obs","sk","ku","pV","rV");
    clear ax ks ku obs pL pV rV sk;
    
    % Practical Salinity (88-21)
    [ax,pL,ks,obs,sk,ku,pV,rV] = L2_ctdHelper(Sp,pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('S_p 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/Sp" + tmpT + ".png");
    save("output\L2\ctd\Sp.mat","pL","ks","obs","sk","ku","pV","rV");
    clear ax ks ku obs pL pV rV sk;
    
    % O2 (88-21)
    [ax,pL,ks,obs,sk,ku,pV,rV] = L2_ctdHelper(o2,pIn,maxMld,dcm,4,"ad",[-60 60]);
    sgtitle('O_2 88-21: L2'+tmpT);
    exportgraphics(ax,"figures/L2/ctd/log/o2" + tmpT + ".png");
    save("output\L2\ctd\o2.mat","pL","ks","obs","sk","ku","pV","rV");
    clear ax ks ku obs pL pV rV sk;

end
%% Unused: NO3-.

% %% No3-: 88-21
% 
% NO3in = nan(329,lowerP,31);
% no3 = nan(lowerP,329);
% 
% for i = cRm
%     tmp = ctdData(i).n(1:lowerP,:);
%     if length(tmp) > 3
%         for j = 1:length(tmp(1,:))
%             NO3in(i,:,j) = tmp(:,j);
%         end
%     end
% end
% 
% for i = 1:329
%     no3(:,i) = mean(squeeze(NO3in(i,:,:)),2,"omitnan");
% end
% 
% [ax,pL,ks,obs,sk,ku,pV,rV] = L2_ctdHelper(no3,pIn,maxMld,dcm);
% sgtitle('NO_3^{-} 88-21: L2');
% exportgraphics(ax,"figures/L2/ctd/log/notUsed/no3" + tmpT + ".png");
% save("output\L2\ctd\no3.mat","pL","ks","obs","sk","ku","pV","rV");
% 
% % ax2 = figure;
% % plot(no3,pL,LineStyle=":",Color=[0.6 0.6 0.6]);
% % set(gca,"YDir","reverse");
% % xlabel('NO_3^{-} [mmol m^{-3}]'); ylabel('P [dbar]'); title('NO_3','P = 0 dbar => DCM');
% % exportgraphics(ax2,'figures/L2/no3.png');
% 
% clearvars -except maxMld dcm lowerP pIn chla T Sp o2 no3 ctdData cRm tmpT;