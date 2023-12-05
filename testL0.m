clear; clc; close all;
addpath("baroneRoutines\");
% addpath("data\L0\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

% Load 
tmp = importdata('data/L0/hplcChla_88-21_200.txt');
pIn = tmp.data(:,4);
chla = tmp.data(:,5);
idIn = tmp.data(:,1);

n = length(pIn);
%% Bin
pB = discretize(pIn,0:10:200);
pX = nan(n,1);

for i = 1:n
    pX(i) = pB(i)*10 - 5;
end

%% K-S
threshold = 50;
n2 = 20;
ks = nan(5,n2);
sk = nan(1,n2); ku = nan(1,n2);
rV = nan(10,n2); pV = nan(10,20);
obs = nan(1,n2);

% 4. Calculate KS p-value, skewness, kurtosis
for i = 1:n2
    % find concentration X_i at binned pressure i
    X_i = chla(pB==i);
    % apply KS test to X_i
    % change limit below to >3 to fix error with picoeu -> may change other
    % results
    if length(X_i) > 3
        [~,ks(:,i),~,~,~,~] = statsplot2(X_i,'noplot');
        [rV(:,i),pV(:,i)] = bbvuong(X_i);
        sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

%%

ks(:,obs<threshold) = nan;

pXX = 5:10:195;
figure
plot(ks(2,:),pXX,'+--','Color','#1f78b4');
set(gca,'YDir','reverse');
ylabel('Pressure [dbar]'); xlabel('p-value'); title('K-S: Lognormal (L0)');