clc;close all;clear;
addpath("baroneRoutines\");

%% CASE 1: Temperature at De Bilt
% 
% tmp = importdata('data/metEx/sum22.txt');
% 
% t = tmp.data(:,2);
% TG = 0.1*tmp.data(:,3);
% % sTG = randsample(TG,250);
% 
% %%
% figure;
% histogram(TG,20);
% 
% %%
% pN = mle(TG,distribution="Normal");
% pG = mle(TG,distribution="Gamma");
% 
% x = (TG-pN(1))/pN(2);
% 
% xCdf = linspace(min(x)-2*std(x),max(x)+2*std(x),2000);
% yCdfN = cdf('norm',xCdf,pN(1),pN(2));
% yCdfG = cdf('logn',xCdf,pG(1),pG(2));
% 
% [h,p,D,cv] = kstest(x);
% [hl,pl,Dl,cvl] = kstest(x,"CDF",[xCdf' yCdfG']);

%% NOTE
% Can't test versus lognormal because of negative temperatures => look at
% wind speed instead to compare e.g. weibull and lognormal.

%% CASE 2: Wind Speed Reanalysis in Hawai'i

finfo = ncinfo("data\metEx\HI_1960_JJA_10mws.nc");
% test = ncread("data\metEx\HI_1960_JJA_10mws.nc",'si10');

%% Load Data

lon = ncread("data\metEx\HI_1960_JJA_10mws.nc","longitude"); 
lat = ncread("data\metEx\HI_1960_JJA_10mws.nc","latitude");
t = ncread("data\metEx\HI_1960_JJA_10mws.nc","time");
ws = ncread("data\metEx\HI_1960_JJA_10mws.nc","si10");  % wind speed (m/s)

wsH0 = squeeze(ws(17,12,:));

nGw = nan(1000,1);
nGl = nan(1000,1);
for i = 1:1000
    wsH = randsample(wsH0,50);
    
    %% K-S
    
    % MLE parameters
    pW = mle(wsH,distribution="Weibull");
    pL = mle(wsH,distribution="Lognormal");
    % pG = mle(wsH,distribution="Gamma");
    pN = mle(wsH,distribution="Normal");
    
    % CDF
    xCdf = linspace(min(wsH)-2*std(wsH),max(wsH)+2*std(wsH),2000);
    yCdfW = cdf('Weibull',xCdf,pW(1),pW(2));
    yCdfL = cdf('Lognormal',xCdf,pL(1),pL(2));
    % yCdfG = cdf('Gamma',xCdf,pG(1),pG(2));
    yCdfN = cdf("Normal",xCdf,pN(1),pN(2));

    % K-S Test
    [hL,pL,DL,cvL] = kstest(wsH,"CDF",[xCdf' yCdfL']);
    [hW,pW,DW,cvW] = kstest(wsH,"CDF",[xCdf' yCdfW']);
    % [hG,pG,DG,cvG] = kstest(wsH,"CDF",[xCdf' yCdfG']);
    [hN,pN,DN,cvN] = kstest(wsH,"CDF",[xCdf' yCdfN']);


%     %% A-D
%     
%     [hLA,pLA,adstatL,cvLA] = adtest(wsH,"Distribution","logn");
%     [hWA,pWA,adstatW,cvWA] = adtest(wsH,"Distribution","weibull");
%     [hNA,pNA,adstatNA,cvNA] = adtest(wsH,"Distribution","norm");

    %% Vuong Ratio Test
    
    [rV,pV] = bbvuong(wsH);
    
    if pN > pW
        nGw(i) = 1;
    end
    if pN > pL
        nGl(i) = 1;
    end
    disp(i);
end

NW = length(find(nGw(nGw==1)))/1000;
NL = length(find(nGl(nGl==1)))/1000;

%% Compare Results

figure
subplot(2,3,1)
histogram(wsH,12);
xlabel('wind speed (m/s)');
subplot(2,3,2)
normplot(wsH);
subplot(2,3,3)
wblplot(wsH);
subplot(2,3,[4 6])
h = cdfplot(wsH);
h.DisplayName = 'Data';
hold on
% plot(xCdf,yCdfL,'DisplayName','Lognormal');
plot(xCdf,yCdfW,'DisplayName','Weibull');
% plot(xCdf,yCdfG,'DisplayName','Gamma');
plot(xCdf,yCdfN,DisplayName='Normal');
hold off
legend();
sgtitle('Wind Speed in Honolulu: 1960 JJA Reanalysis');
