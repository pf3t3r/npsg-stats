close all;clc;clear;
addpath 'C:\Users\pfarrell\AppData\Roaming\MathWorks\MATLAB Add-Ons\Functions\Shapiro-Wilk and Shapiro-Francia normality tests'



%%
% figure;
% histogram(tmp,25);

%%
r = 300;
s = 250;

for i = 1:r
    tmp = randn(s,1);
    [h(i),p(i)] = kstest(tmp);
    [h1(i),p1(i)] = lillietest(tmp);
    [h2(i),p2(i)] = adtest(tmp);
    [h3(i),p3(i)] = swtest(tmp);
    [h4(i),p4(i)] = chi2gof(tmp);
    [h5(i),p5(i)] = jbtest(tmp);
end

% clear h h1 h2 h3 h4 h5 tmp;
H = [h h1 h2 h3 h4 h5];

% Percentage of incorrect results
% i.e. where a test reports that the distribution is NOT normal
X = length(find(H(H==1)))/length(H);