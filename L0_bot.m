% Script to output L0 bottle results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\");
addpath("func\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
% tmpT = "";

% Extra possible test cases.
analyseYearByYear = false;  % analyse effect (if any) of varying start year 
                            % on distributions
crn131 = false;             % analyse 2001-2021 data (to mirror CTD results)
nightAnalysis = false;      % analyse night-time 2001-2021 (to mirror CTD 
                            % results)
logAxes = false;            % output p-values as log values

if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% K-S: Principal Analysis.

tmpT = "";

% chla
tmp = importdata('data/L0/hplcChla_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Chl $a$ (1988-2021)","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% mvchla 
tmp = importdata('data/L0/mvchla_94-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Monovinyl Chl $a$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% dvchla
tmp = importdata('data/L0/dvchla_94-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Divinyl Chl $a$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% chlb
tmp = importdata('data/L0/chlb_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Chl $b$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% chlc123
tmp = importdata('data/L0/chlc123_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Chl $c1 + c2 + c3$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% alpha-caro
tmp = importdata('data/L0/acar_94-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: $\alpha$-carotene","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% but-19
tmp = importdata('data/L0/but19_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: 19' Butanoyloxyfucoxanthin","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% hex-19
tmp = importdata('data/L0/hex19_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: 19' Hexanoyloxyfucoxanthin","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% zeax
tmp = importdata('data/L0/zeax_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Zeaxanthin","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% pc
tmp = importdata('data/L0/pc_89-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Particulate Carbon","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% pn
tmp = importdata('data/L0/pn_89-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Particulate Nitrogen","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% llp
tmp = importdata('data/L0/llp_88-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Low-Level Phosphorus","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% lln
tmp = importdata('data/L0/lln_88-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Low-Level Nitrogen","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% nit + nit
tmp = importdata('data/L0/nitNit_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Nitrate + Nitrite","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% phos
tmp = importdata('data/L0/phos_88-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Phosphate","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% boxy
tmp = importdata('data/L0/boxy_88-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Dissolved Oxygen","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% L-12 PP
tmp = importdata('data/L0/l12_89-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
sgtitle("L0: Primary Production (Light-12)","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

%% A-D: Principal Analysis.

tmpT = "-ad";

% chla
tmp = importdata('data/L0/hplcChla_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "chla" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% mvchla 
tmp = importdata('data/L0/mvchla_94-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% dvchla
tmp = importdata('data/L0/dvchla_94-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% chlb
tmp = importdata('data/L0/chlb_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% chlc123
tmp = importdata('data/L0/chlc123_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% alpha-caro
tmp = importdata('data/L0/acar_94-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% but-19
tmp = importdata('data/L0/but19_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% hex-19
tmp = importdata('data/L0/hex19_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% zeax
tmp = importdata('data/L0/zeax_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% pc
tmp = importdata('data/L0/pc_89-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% pn
tmp = importdata('data/L0/pn_89-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% llp
tmp = importdata('data/L0/llp_88-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% lln
tmp = importdata('data/L0/lln_88-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% nit + nit
tmp = importdata('data/L0/nitNit_88-21_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% phos
tmp = importdata('data/L0/phos_88-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% boxy
tmp = importdata('data/L0/boxy_88-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

% L-12 PP
tmp = importdata('data/L0/l12_89-22_200.txt');
[ax,~,~] = L0_helper(tmp,50,'ad');
sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

%% 2001-2021 chla (starting CRN 131, to match newer fluorometer)
if crn131 == true
    
    tmpT = "";
    
    % K-S
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2001-2021)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21" + tmpT + ".png");
    clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;
    
    % A-D
    tmpT = "-ad";
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2001-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21" + tmpT + ".png");
    clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

else
    disp("Not analysing CRN 131 only data...");
end

%% Chl-a 2001-2022 NIGHT-TIME
if nightAnalysis == true
        
    tmpT = "";
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    
    hms = char(string(tmp.data(:,3)));
    mdy = char(string(tmp.data(:,2)));
    
    % test = "0" + t2(1,1:end-1);
    for i = 1:length(tmp.data)
        if hms(i,end) == " "
            hms(i,:) = "0" + hms(i,1:end-1);
        end
        if mdy(i,end) == " "
            mdy(i,:) = "0" + mdy(i,1:end-1);
        end
    end
    
    Y = double("20" + mdy(:,5:6));
    M = double("" + mdy(:,1:2));
    D = double("" + mdy(:,3:4));
    h = double("" + hms(:,1:2));
    m = double("" + hms(:,3:4));
    s = double("" + hms(:,5:6));
    
    T = datetime(Y,M,D,h,m,s);
    T2 = datetime(Y,M,D);
    T3 = datenum(T);
    
    Lat = 22.75;
    Lon = -158;
    [SunRiseSet,~,~,~,~,~] = suncycle(Lat,Lon,T2);
    for i = 1:length(SunRiseSet)
        tmp = SunRiseSet(i,:) - 10;
        for j = 1:2
            if tmp(j) < 0
                tmp(j) = tmp(j) + 24;
            end
        end
        rs(i,:) = tmp;
    end
    rs2 = hours(rs);
    rs2.Format = 'hh:mm';
    
    sunriseTime = hours(rs2(:,1)');
    sunsetTime = hours(rs2(:,2)');
    
    dayFrac = rem(T3,1);
    castTime = dayFrac*24;
    
    castAtNight = nan(length(T),1);
    
    for i = 1:length(T)
        if castTime(i) < sunriseTime(i) || castTime(i) > sunsetTime(i)
            castAtNight(i) = 1;
        end
    end
    
    nightCastIDs = [];
    
    for i=1:length(T)
        if(~isnan(castAtNight(i)))
            disp(i);
            nightCastIDs = [nightCastIDs i];
        end
    end
    
    tmp = importdata('data/L0/hplcChla_01-22_200.txt').data(nightCastIDs,:);
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2001-2021, NIGHT)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21_night" + tmpT + ".png");
    clearvars -except tmpT analyseYearByYear crn131 nightAnalysis logAxes lp;

else
    disp("Not analysing night-time only data...");
end

%% Year-by-year analysis for chla
% Here we move the start date of the analysis forward in time to see if
% the distribution of data has some dependence on time. It will only be
% checked if "analyseYearByYear" is true.
if analyseYearByYear==true
    
    % K-S
    tmpT = "";
    
    % 2016-2022
    tmp = importdata('data/L0/hplcChla_16-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2016-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_16-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2015-2022
    tmp = importdata('data/L0/hplcChla_15-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2015-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_15-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2013-2022
    tmp = importdata('data/L0/hplcChla_13-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2013-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_13-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2012-2022
    tmp = importdata('data/L0/hplcChla_12-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2012-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_12-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2011-2022
    tmp = importdata('data/L0/hplcChla_11-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2011-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_11-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2011-2022
    tmp = importdata('data/L0/hplcChla_11-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2011-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_11-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2010-2022
    tmp = importdata('data/L0/hplcChla_10-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2010-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_10-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2009-2022
    tmp = importdata('data/L0/hplcChla_09-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2009-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_09-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2008-2022
    tmp = importdata('data/L0/hplcChla_08-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2008-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_08-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2007-2022
    tmp = importdata('data/L0/hplcChla_07-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2007-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_07-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2006-2022
    tmp = importdata('data/L0/hplcChla_06-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2006-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_06-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2005-2022
    tmp = importdata('data/L0/hplcChla_05-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2005-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_05-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2004-2022
    tmp = importdata('data/L0/hplcChla_04-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2004-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_04-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2003-2022
    tmp = importdata('data/L0/hplcChla_03-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2003-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_03-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2002-2022
    tmp = importdata('data/L0/hplcChla_02-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2002-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_02-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2001-2022 done above.
    
    % 2000-2022
    tmp = importdata('data/L0/hplcChla_00-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2000-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_00-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1999-2022
    tmp = importdata('data/L0/hplcChla_99-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1999-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_99-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1998-2022
    tmp = importdata('data/L0/hplcChla_98-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1998-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_98-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1997-2022
    tmp = importdata('data/L0/hplcChla_97-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1997-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_97-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1996-2022
    tmp = importdata('data/L0/hplcChla_96-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1996-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_96-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1995-2022
    tmp = importdata('data/L0/hplcChla_95-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1995-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_95-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1994-2022
    tmp = importdata('data/L0/hplcChla_94-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1994-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_94-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1993-2022
    tmp = importdata('data/L0/hplcChla_93-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1993-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_93-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1992-2022
    tmp = importdata('data/L0/hplcChla_92-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1992-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_92-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1991-2022
    tmp = importdata('data/L0/hplcChla_91-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1991-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_91-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1990-2022
    tmp = importdata('data/L0/hplcChla_90-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1990-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_90-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1989-2022
    tmp = importdata('data/L0/hplcChla_89-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1989-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_89-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    
    %% chla (A-D): go back year by year.
    tmpT = "-ad";
    
    % 2016-2022
    tmp = importdata('data/L0/hplcChla_16-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2016-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_16-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2015-2022
    tmp = importdata('data/L0/hplcChla_15-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2015-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_15-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2014-2022
    tmp = importdata('data/L0/hplcChla_14-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2014-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_14-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2013-2022
    tmp = importdata('data/L0/hplcChla_13-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2013-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_13-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2012-2022
    tmp = importdata('data/L0/hplcChla_12-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2012-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_12-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2011-2022
    tmp = importdata('data/L0/hplcChla_11-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2011-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_11-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2010-2022
    tmp = importdata('data/L0/hplcChla_10-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2010-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_10-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2009-2022
    tmp = importdata('data/L0/hplcChla_09-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2009-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_09-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2008-2022
    tmp = importdata('data/L0/hplcChla_08-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2008-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_08-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2007-2022
    tmp = importdata('data/L0/hplcChla_07-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2007-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_07-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2006-2022
    tmp = importdata('data/L0/hplcChla_06-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2006-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_06-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2005-2022
    tmp = importdata('data/L0/hplcChla_05-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2005-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_05-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2004-2022
    tmp = importdata('data/L0/hplcChla_04-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2004-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_04-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2003-2022
    tmp = importdata('data/L0/hplcChla_03-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2003-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_03-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2002-2022
    tmp = importdata('data/L0/hplcChla_02-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2002-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_02-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 2001-2022 done above.
    
    % 2000-2022
    tmp = importdata('data/L0/hplcChla_00-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2000-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_00-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1999-2022
    tmp = importdata('data/L0/hplcChla_99-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1999-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_99-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1998-2022
    tmp = importdata('data/L0/hplcChla_98-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1998-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_98-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1997-2022
    tmp = importdata('data/L0/hplcChla_97-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1997-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_97-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1996-2022
    tmp = importdata('data/L0/hplcChla_96-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1996-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_96-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1995-2022
    tmp = importdata('data/L0/hplcChla_95-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1995-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_95-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1994-2022
    tmp = importdata('data/L0/hplcChla_94-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1994-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_94-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1993-2022
    tmp = importdata('data/L0/hplcChla_93-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1993-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_93-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1992-2022
    tmp = importdata('data/L0/hplcChla_92-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1992-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_92-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1991-2022
    tmp = importdata('data/L0/hplcChla_91-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1991-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_91-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1990-2022
    tmp = importdata('data/L0/hplcChla_90-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1990-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_90-22" + tmpT + ".png");
    clearvars -except tmpT;
    
    % 1989-2022
    tmp = importdata('data/L0/hplcChla_89-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1989-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_89-22" + tmpT + ".png");
    clearvars -except tmpT;

else
    disp("Not analysing year by year...");
end