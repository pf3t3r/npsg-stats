close all;clc;clear;

tmp = importdata("data\L0\t_88-22-200.txt").data;

crn = tmp(:,1);
day = tmp(:,2);
p = tmp(:,3);
T = tmp(:,4);

%% Overall
figure;
scatter(T,p,13,"green",'.');
set(gca,"YDir","reverse");
ylabel("Pressure [dbar]"); xlabel("T (C)");

%% Cruise 1
a = 1; b=101;
figure;
scatter(T(a:b),p(a:b),10,"red",'.');
set(gca,"YDir","reverse");
ylabel("Pressure [dbar]"); xlabel("T (C)");

%% Cruise 100
a = 9798; b=9898;
figure;
scatter(T(a:b),p(a:b),10,"red",'.');
set(gca,"YDir","reverse");
ylabel("Pressure [dbar]"); xlabel("T (C)");

%%

% figure;
% normplot(T);
T2 = T(T>0);

figure;
probplot("lognormal",T2);