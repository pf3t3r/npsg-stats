close all; clc; clear;
addpath("func\");

% I am doing this wrong. I need to consider the depths. o_o

%% Load Sample Data
chl_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,5);
cmo = importdata('data/chlaMonovinyl_88-21.txt').data(:,5);
car = importdata('data/parC_89-20.txt').data(:,5);

% Remove zeros or negative values
chl_hplc(chl_hplc<=0) = nan;
cmo(cmo<=0) = nan;
car(car<=0) = nan;

%% Visualise distributions for different variables

funcVisualiseDistributions(chl_hplc,'Chl a (HPLC)');
funcVisualiseDistributions(cmo,'Monovinyl Chl a (HPLC)');
% funcVisualiseDistributions(car,'Particulate Carbon');

%%
p_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,4);
% chl_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,5);
id_hplc = importdata('data/HPLC_chla_88-21.txt').data(:,1);
[pb5_hplc,pb10_hplc,chlOut_hplc,n5_hplc,n10_hplc] = cleanAndBin(p_hplc,chl_hplc,id_hplc);

ks = nan(5,20); obs = nan(20,1); n = 20; depth = 5:10:200;
for i = 10:10
    X_i = chl_hplc(pb10_hplc==i);
    funcVisualiseDistributions(X_i,'Chl a (HPLC): 100 dbar');
end


%% Reference: Examples of Distribution Types

% biSa = makedist("BirnbaumSaunders");
% % xBisa = 0:20;
% % yBisa = pdf(biSa,xBisa);
% burr = makedist("Burr");
% expo = makedist("Exponential");
% exVa = makedist("ExtremeValue");
% gamm = makedist("Gamma");
% geEx = makedist("GeneralizedExtremeValue");
% gePa = makedist("GeneralizedPareto");
% haNo = makedist("HalfNormal");
% inGa = makedist("InverseGaussian");
% logi = makedist("Logistic");
% logl = makedist("Loglogistic");
% logn = makedist("Lognormal");
% naka = makedist("Nakagami");
% norm = makedist("Normal");
% pois = makedist("Poisson");
% rayl = makedist("Rayleigh");
% rici = makedist("Rician");
% tLoc = makedist("tLocationScale");
% unif = makedist("Uniform");
% weib = makedist("Weibull");
% 
% figure;
% % plot(xBisa,yBisa);
% plot(biSa);
% hold on
% plot(burr); plot(expo); plot(exVa); plot(gamm); plot(geEx); plot(gePa);
% plot(haNo); plot(inGa); plot(logi); plot(logl); plot(naka); plot(norm);
% plot(pois); plot(rayl); plot(rici); plot(tLoc); plot(unif); plot(weib);
% hold off
% legend();
% xlim([0 10]);
