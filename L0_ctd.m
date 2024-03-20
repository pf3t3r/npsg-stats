% Script to output L0 ctd results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

% %% CTD casts where hplc chl-a samples taken
% 
% % CRN 131 - 339?
% 
% % !!!! what am I doing here!!!
% data = importdata("data\L0\hplcChla_01-22_200.txt").data;
% botid = char(string(data(:,1)));
% crn = str2num(botid(:,1:3));
% cast = str2num(botid(:,6:8));
% 
% crnCast = [crn cast];
% 
% crnCastCombo = unique(crnCast,"rows");
% id = 1:1:200;
% 
% crnCastCombo = [id' crnCastCombo];
% 
% crn2 = unique(crn);
% 
% fuck = load("datafiles\ctd_iso_ALL.mat").ctd;
% j = 0;
% 
% % CRN 181
% % f181 = mean(fuck(181).f(1:101,[8 15]),2);
% 
% for i = 1:191
%     disp(crnCastCombo(i,:));
%     tmp = crnCastCombo(i,:);
%     f(:,i) = fuck(tmp(2)).f(1:101,tmp(3));
%     time(i) = fuck(tmp(2)).decimal_hour(tmp(3));
% end
% 
% time = time - 10;
% for i = 1:length(time)
%     if time(i) < 0
%         time(i) = time(i) + 24;
%     end
% end
%% CTD data

% This data is mean of all casts. Downloaded via FTP and averaged in other
% file. It SHOULD be equivalent to the CTD chloropigment we find on the
% HOT-DOGS system.

%% K-S
% 2001-2021 (crn 131-)
tmpT = "";
chla = load("output\CTD\chla.mat").meanEpN(1:101,131:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2002-2021 (crn 134-)
tmpT = "02";
chla = load("output\CTD\chla.mat").meanEpN(1:101,134:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2003-2021 (crn 144-)
tmpT = "03";
chla = load("output\CTD\chla.mat").meanEpN(1:101,144:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2004-2021 (crn 155-)
tmpT = "04";
chla = load("output\CTD\chla.mat").meanEpN(1:101,155:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2005-2021 (crn 167-)
tmpT = "05";
chla = load("output\CTD\chla.mat").meanEpN(1:101,167:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2006-2021 (crn 177-)
tmpT = "06";
chla = load("output\CTD\chla.mat").meanEpN(1:101,177:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2007-2021 (crn 189-)
tmpT = "07";
chla = load("output\CTD\chla.mat").meanEpN(1:101,189:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2008-2021 (crn 199-)
tmpT = "08";
chla = load("output\CTD\chla.mat").meanEpN(1:101,199:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2009-2021 (crn 208-)
tmpT = "09";
chla = load("output\CTD\chla.mat").meanEpN(1:101,208:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2010-2021 (crn 219-)
tmpT = "10";
chla = load("output\CTD\chla.mat").meanEpN(1:101,219:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2011-2021 (crn 228-)
tmpT = "11";
chla = load("output\CTD\chla.mat").meanEpN(1:101,228:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2012-2021 (crn 239-)
tmpT = "12";
chla = load("output\CTD\chla.mat").meanEpN(1:101,239:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2013-2021 (crn 249-)
tmpT = "13";
chla = load("output\CTD\chla.mat").meanEpN(1:101,249:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2014-2021 (crn 259-)
tmpT = "14";
chla = load("output\CTD\chla.mat").meanEpN(1:101,259:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2015-2021 (crn 269-)
tmpT = "15";
chla = load("output\CTD\chla.mat").meanEpN(1:101,269:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2016-2021 (crn 280-)
tmpT = "16";
chla = load("output\CTD\chla.mat").meanEpN(1:101,280:329);
ax = L0_ctdHelper(chla);
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

%% A-D
% 2001-2021 (crn 131-)
tmpT = "ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,131:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2002-2021 (crn 134-)
tmpT = "02-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,134:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2003-2021 (crn 144-)
tmpT = "03-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,144:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2004-2021 (crn 155-)
tmpT = "04-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,155:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2005-2021 (crn 167-)
tmpT = "05-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,167:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2006-2021 (crn 177-)
tmpT = "06-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,177:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2007-2021 (crn 189-)
tmpT = "07-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,189:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2008-2021 (crn 199-)
tmpT = "08-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,199:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2009-2021 (crn 208-)
tmpT = "09-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,208:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2010-2021 (crn 219-)
tmpT = "10-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,219:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2011-2021 (crn 228-)
tmpT = "11-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,228:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2012-2021 (crn 239-)
tmpT = "12-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,239:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2013-2021 (crn 249-)
tmpT = "13-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,249:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2014-2021 (crn 259-)
tmpT = "14-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,259:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2015-2021 (crn 269-)
tmpT = "15-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,269:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;

% 2016-2021 (crn 280-)
tmpT = "16-ad";
chla = load("output\CTD\chla.mat").meanEpN(1:101,280:329);
ax = L0_ctdHelper(chla,"ad");
sgtitle("chl-a " + tmpT);
exportgraphics(ax,"figures/L0/ctd/chla" + tmpT + ".png"); clear;
%%
% close all;
% 
% tmpAlpha = importdata('data/L0/alphaCaro_50-60.txt');
% tmpBut19 = importdata('data\L0\but19_50-60.txt');
% tmpHex19 = importdata('data\L0\hex19_50-60.txt');
% 
% alphaCar = tmpAlpha.data(:,5);
% but19 = tmpBut19.data(:,5);
% hex19 = tmpHex19.data(:,5);
% 
% %%
% figure; histogram(alphaCar,max(alphaCar)); title('alpha-carotene');
% 
% figure; histogram(but19,max(but19)); title('But-19');
% 
% figure; histogram(hex19,max(hex19)); title('Hex-19');
% 
% %%
% [mleA,ksA] = statsplot2(alphaCar);
% [~,ksB] = statsplot2(but19);
% [~,ksH] = statsplot2(hex19);
% 
% 
% pd = makedist("Gamma","a",mleA(4,1),"b",mleA(4,2));
% 
% [~,pAdA] = adtest(alphaCar,"Distribution",pd);
% [~,pAdB] = adtest(but19,"Distribution",pd);
% [~,pAdH] = adtest(hex19,"Distribution",pd);
% 
% [rA,pA] = bbvuong(alphaCar);
% [rB,pB] = bbvuong(but19);
% [rH,pH] = bbvuong(hex19);
% 
% %%
% dAcar = load("output\L1\acar.mat");
% pOutA = dAcar.pOut;
% cOutA = dAcar.cOut;
% 
% dHex19 = load("output\L1\hex19.mat");
% pOutH = dHex19.pOut;
% cOutH = dHex19.cOut;
% 
% dBut19 = load("output\L1\but19.mat");
% pOutB = dBut19.pOut;
% cOutB = dBut19.cOut;
% %%
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 36 10]);
% 
% figure;
% subplot(1,3,1)
% histogram(alphaCar); title("alpha-caro");
% subplot(1,3,2)
% histogram(but19); title("but-19");
% subplot(1,3,3)
% histogram(hex19); title("hex-19");
% sgtitle('50-60 dbar bin: L0');
% 
% %%
% figure;
% subplot(1,3,1)
% histogram(cOutA(pOutA==6)); title("alpha-caro");
% subplot(1,3,2)
% histogram(cOutB(pOutB==6)); title("but-19");
% subplot(1,3,3)
% histogram(cOutH(pOutH==6)); title("hex-19");
% sgtitle('50-60 dbar bin: L1');
% 
% %% remove outlier in hex-19 L0
% 
% hex19o = hex19;
% hex19o(40) = [];
% % show
% figure
% plot(hex19);
% hold on
% plot(hex19o);
% hold off
% 
% [~,ksHo] = statsplot2(hex19o);
% [rHo,pHo] = bbvuong(hex19o);