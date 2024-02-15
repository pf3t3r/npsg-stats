clear; clc; close all;
addpath("baroneRoutines\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

tmpT = "";
%%
tmp = importdata('data/L0/hplcChla_88-21_200.txt');
[ax,~,~] = L0_helper(tmp);
% sgtitle("L0: Chl $a$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/chla" + tmpT + ".png");
clearvars -except tmpT;

%% Chl-a: 2001-2021 (from CRN 131, to match newer fluorometer)
tmp = importdata('data/L0/hplcChla_01-22_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Chl $a$ (2001-2021)","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/chla_01-21" + tmpT + ".png");
clearvars -except tmpT;

%% Chl-a: 2001-2022, NIGHT-TIME

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

tmp = importdata('data/L0/pc_89-22_200.txt').data(nightCastIDs,:);
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Chl $a$ (2001-2021, NIGHT)","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/chla_01-21_night" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/pc_89-22_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Particulate Carbon","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/pc" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/mvchla_94-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Monovinyl Chl $a$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/mvchla" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/dvchla_94-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Divinyl Chl $a$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/dvchla" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/chlb_88-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Chl $b$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/chlb" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/chlc123_88-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Chl $c1 + c2 + c3$","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/chlc123" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/acar_94-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: $\alpha$-carotene","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/acar" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/but19_88-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: 19' Butanoyloxyfucoxanthin","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/but19" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/hex19_88-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: 19' Hexanoyloxyfucoxanthin","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/hex19" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/zeax_88-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Zeaxanthin","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/zeax" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/pn_89-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Particulate Nitrogen","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/pn" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/boxy_88-22_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Dissolved Oxygen","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/boxy" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/lln_88-22_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Low-Level Nitrogen","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/lln" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/nitNit_88-21_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Nitrate + Nitrite","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/nit2" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/llp_88-22_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Low-Level Phosphorus","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/llp" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/phos_88-22_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Phosphate","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/phos" + tmpT + ".png");
clearvars -except tmpT;
%%
tmp = importdata('data/L0/l12_89-22_200.txt');
[ax,~,~] = L0_helper(tmp);
sgtitle("L0: Primary Production (Light-12)","Interpreter","latex");
exportgraphics(ax,"figures/L0/bot/l12" + tmpT + ".png");
clearvars -except tmpT;