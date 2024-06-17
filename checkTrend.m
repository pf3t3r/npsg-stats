% Check the trend in chl-a in 1988-2021

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 28 15]);

% Load maximum mixed-layer depth 'MLD' and cruise-averaged deep chlorophyll
% maximum depth 'DCM'.
mld = load("mldVals.mat").maxMld; % single maximum per cruise
dcm = load("output/dcm.mat").dcm; % pDcm + sigmaDcm (all casts, crn 1-329)

%  K-S
tmpT = "";

% Chlorophyll a (88-21)
tmp = importdata("data/L2/hplcChla_88-21_200.txt");
MDY = tmp.data(:,2);
MDY(MDY == -9) = 999944;
MDY = num2str(MDY);
hms = tmp.data(:,3);
hms(hms == -9) = 999999;
hms = num2str(hms);

% for i = 1:length(MDY)
%     if ~isnan(MDY(i))
%         disp(i);
%         M = str2num(MDY(:,1:2)); D = str2num(MDY(:,3:4)); Y = str2num(MDY(:,5:6));
%     end
% end

for i = 1:length(MDY)
    MDY(i,:) = strrep(MDY(i,:),' ','0');
end
M = str2num(MDY(:,1:2)); D = str2num(MDY(:,3:4)); Y = str2num(MDY(:,5:6));

for i = 1:length(MDY)
    if Y(i) > 87
        tmp = strcat('19',num2str(Y(i)));
        Y(i) = str2num(tmp);
    elseif Y(i) < 10
        tmp = strcat('200',num2str(Y(i)));
        Y(i) = str2num(tmp);
    elseif Y(i) < 87
        tmp = strcat('20',num2str(Y(i)));
        Y(i) = str2num(tmp);
    end
end

for i = 1:length(hms)
    hms(i,:) = strrep(hms(i,:),' ','0');
end
h = str2num(hms(:,1:2)); m = str2num(hms(:,3:4)); s = str2num(hms(:,5:6));
t = datetime(Y,M,D,h,m,s);

tmp = importdata("data/L2/hplcChla_88-21_200.txt");
[ax,X,p,bA] = L2_helper(tmp,mld,dcm,50,4,"ks",[-60 60],[7 19]);
sgtitle("[Chl a] 88-21: L2");
exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");

binnedPressure = round(bA(:,5),-1);
figure;
scatter(X,binnedPressure,'.');
set(gca,"YDir","reverse");

% 0 line
p0 = binnedPressure(binnedPressure==0);
x0 = X(binnedPressure==0);

figure;
plot(x0);

% LT = trenddecomp(x0);
% figure;
% plot()