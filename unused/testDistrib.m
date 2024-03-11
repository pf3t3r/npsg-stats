clc; close all; clear;

% For this example we load data for phosphate from between 70 and 80 dbar.
% At this depth we expect the normal and weibull distributions to be the
% most important (of the initial five sample distributions).

% ISSUE: The p-values are much smaller than expected.

%% Load Phosphate (70-80 dbar, 1988-2021)
phoId = num2str(importdata('data/dist/pho70m.txt').data(:,1));
phoP = importdata('data/dist/pho70m.txt').data(:,4);
pho = importdata('data/dist/pho70m.txt').data(:,5);

% Remove zero values as this causes problems for some distributions,
% notably Birnbaum-Saunders.
pho(pho==0) = NaN;

sil = importdata('data/dist/sil70.txt').data(:,5);
%%
figure; histogram(pho,"BinEdges",linspace(0.01,0.16,16));

% Low precision. Phosphorus takes only one of 16 possible values.
% min = 0.01, max = 0.16, mean = 0.06.

%%
figure; histogram(sil);
uniquePho = length(unique(pho));
uniqueSil = length(unique(sil));
%% Find MLE

% Distributions used:
% Normal
% Distributions not used / not appropriate:
% Bernouilli, 
X = sil;

MLEp = nan(6,2); KSp = nan(6,1);

MLEp(1,:) = mle(X);
% MLEp(2,:) = mle(X,"distribution","Beta");
MLEp(3,1) = mle(X,'distribution','Exponential');
MLEp(4,:) = mle(X,'distribution','Extreme Value');
MLEp(5,:) = mle(X,distribution="Gamma");
MLEp(6,:) = mle(X,distribution="Lognormal");
% MLEp(3,:) = mle(pho,'distribution','Binomial','NTrials',100); % no valid rows?
%MLEp(3,:) = mle(pho,'distribution','birnbaumsaunders');
% MLEp(5,:) = mle(pho,'distribution','Burr'); % BURR needs three outputs,
% add later
%MLEp(4,:) = mle(pho,'distribution','Discrete Uniform');


%% Find KS
xCdf = linspace(min(pho)-2*std(pho,'omitnan'),max(pho)+2*std(pho,'omitnan'),2000);
yCdfNorm = cdf('norm',xCdf,MLEp(1,1),MLEp(1,2));
% yCdfBeta = cdf("Beta",xCdf,MLEp(2,1),MLEp(2,2));
yCdfExp = cdf("Exponential",xCdf,MLEp(3,1));
yCdfExt = cdf("Extreme Value",xCdf,MLEp(4,1),MLEp(4,2));
yCdfGam = cdf("Gamma",xCdf,MLEp(5,1),MLEp(5,2));
yCdfLog = cdf("Lognormal",xCdf,MLEp(6,1),MLEp(6,2));
% yCdfBino = cdf("Binomial",xCdf,MLEp(3,1),MLEp(3,2));
%yCdfBirn = cdf("BirnbaumSaunders",xCdf,MLEp(4,1),MLEp(4,2));
%yCdfDisc = cdf("Discrete Uniform",xCdf,MLEp(5,1),MLEp(5,2));

[~,KSp(1)] = kstest(pho,[xCdf' yCdfNorm']);
if sum(pho<=0) == 0
    %[~,KSp(2)] = kstest(pho,[xCdf' yCdfBeta']);
    [~,KSp(3)] = kstest(pho,[xCdf' yCdfExp']);
    [~,KSp(4)] = kstest(pho,[xCdf' yCdfExt']);
    [~,KSp(5)] = kstest(pho,[xCdf' yCdfGam']);
    [~,KSp(6)] = kstest(pho,[xCdf' yCdfGam']);
end

%%
figure;
scatter(KSp,70);
legend();

% Best distributions here: Beta, Normal, Lognormal...
% But p-values way too low.