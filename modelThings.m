clear; clc; close all; addpath("func\");

% Initial Conditions
x0 = 1e-3;          % initial concentration
t0 = 0;             % start time
tF = 100;          % finish time
h = tF/10000;       % step size
time = t0:h:tF;     % time vector

%% Case 1: "Campbell"
% d[X]/dt = -k[X]

k = 1;              % net growth rate
F = @(t,x) -k*x;    % Eqn 1 in Andersson (2021); main equation in Campbell (1995)

C1 = modelEquation(F,t0,h,tF,x0);

%% Case 2: "Andersson". Deterministic + stochastic component to growth.
% Andersson (2021) proposes splitting the net growth rate k into a
% deterministic, constant mean component mu_k, and a stochastic fluctuation
% described by sigma_k*eta(t). Sigma_k is the magnitude of the stochastic
% fluctuation and eta(t) describes the time-dependency of the random
% fluctuations.
% d[X]/dt = -(mu_k + sigma_k*eta(t))[X]

mu_k = 0.01;           % deterministic component of growth k
sigma_k = 1;      % magnitude of the stochastic fluctuation

% default tolerances: RelTol = 1e-3, AbsTol = 1e-6.
[~,C2] = ode45(@(t,x) andersson2(x,mu_k,sigma_k),time,x0);

%% Case 4. Add in the mixing as a constant.
% d[X]/dt = -(mu_k + sigma_k*eta(t))[X] + M.

M = 1e-4;           % mixing term

[t3,C3] = ode45(@(t,x) andersson2(x,mu_k,sigma_k,M),time,x0);

%% Case 5. Mixing varies with time.
% d[X]/dt = -(mu_k + sigma_k*eta(t))[X] + M(t).

[t4,C4] = ode45(@(t,x) andersson2(x,mu_k,sigma_k,M,true),time,x0);

%% Compare Cases 1 - 3.
% Compare the three cases visually.

figure
subplot(4,2,[1 3 5 7])
plot(time,C1,DisplayName="Eq. 1"); hold on
plot(time,C2,DisplayName="Eq. 2");
plot(time,C3,DisplayName="Eq. 2 + constant mixing");
plot(time,C4,DisplayName="Eq. 2 + time-varying mixing");
hold off
grid on
xlabel("time"); ylabel("concentration"); legend();
subplot(4,2,[2 4 6 8])
histogram(C1,50,DisplayName="Eq. 1")
hold on
histogram(C2,50,DisplayName="Eq. 2");
histogram(C3,50,DisplayName="Eq. 2 + constant mixing");
histogram(C4,50,DisplayName="Eq. 2 + time-varying mixing");
hold off
legend();
sgtitle("k = " + k + ", \mu = "+mu_k+", \sigma_k = "+sigma_k+"");
