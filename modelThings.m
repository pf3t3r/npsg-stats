clear; clc; addpath("func\"); close all; 

% Initial Conditions
x0 = 1;             % initial concentration
t0 = 0;             % start time
tF = 10*1000;          % finish time
% h = tF/100000;      % step size
h = 0.01;
time = t0:h:tF;     % time vector

% Cases 1-2 will be run by default. Other cases may be run or not.
runCase3 = true;
runCase4 = true;
checkEtas = false;

%% Case 1: "Campbell"
% d[X]/dt = -k[X]
% For the purposes of this model we rename the variable 'k' from Campbell
% (1995) to mu_k in order to be consistent with the equation from Andersson
% (2021).

% mu_k = -10/tF;               % net growth rate / deterministic component of growth
mu_k = -0.00001;
% mu_k = 0;
F = @(t,x) -mu_k*x;         % Eqn 1 in Andersson (2021); main equation in Campbell (1995)

C1 = modelEquation(F,t0,h,tF,x0);

%% Case 2: "Andersson". Deterministic + stochastic component to growth.
% Andersson (2021) proposes splitting the net growth rate k into a
% deterministic, constant mean component mu_k, and a stochastic fluctuation
% described by sigma_k*eta(t). Sigma_k is the magnitude of the stochastic
% fluctuation and eta(t) describes the time-dependency of the random
% fluctuations.
% d[X]/dt = -(mu_k + sigma_k*eta(t))[X]

sigma_k = 10*abs(mu_k);             % magnitude of the stochastic fluctuation
% sigma_k = 5e-3;

opt = odeset(Events=@stopIntegration);  % stop integration at 0.001*x0

C2 = nan(length(time),5); t2 = nan(length(time),5);
for i = 1:5
    [tmpT,tmp] = ode45(@(t,x) andersson2(x,mu_k,sigma_k),time,x0,opt);
    C2(1:length(tmp),i) = tmp;
    t2(1:length(tmpT),i) = tmpT;
end

if checkEtas == true
    for i = 1:5
        for k = 1:numel(t2(:,1))
            [~,etas(k,i)] = andersson2(C2(k,i),mu_k,sigma_k);
        end
    end
end

%% Case 3. Add in the mixing as a constant.
% d[X]/dt = -(mu_k + sigma_k*eta(t))[X] + M.

C3 = nan(length(time),5);
if runCase3 == true
    xTarget = 0.05;
    %mix = xTarget*mu_k;           % mixing term
    mix = -1e-5;

    for i = 1:5
        [~,tmp] = ode45(@(t,x) andersson2(x,mu_k,sigma_k,mix),time,x0,opt);
        C3(1:length(tmp),i) = tmp;
    end
end
%% Case 4. Mixing varies with time.
% d[X]/dt = -(mu_k + sigma_k*eta(t))[X] + M(t).

C4 = nan(length(time),5);
if runCase4 == true
    for i = 1:5
        [~,tmp] = ode45(@(t,x) andersson2(x,mu_k,sigma_k,mix,true),time,x0,opt);
        C4(1:length(tmp),i) = tmp;
    end
end
%% Compare Cases 1 - 4 visually.

% Plot also theoretical skewness/kurtosis for a lognormal distribution
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

figure
subplot(1,3,1)
plot(time,C1,LineStyle="-.",LineWidth=1,Color="#d95f02"); hold on
plot(time,C2,Color='#7570b3');
if runCase3 == true
    plot(time,C3,Color='#e7298a');
end
if runCase4 == true
    plot(time,C4,Color='#66a61e');
end
hold off; grid on
xlabel("time"); ylabel("concentration");
legend("Eq. 1","Eq.2","","","","","Eq. 2 + constant mixing","","","","","Eq. 2 + time-varying mixing","","","","");

subplot(1,3,2)
histogram(C1,50,DisplayName="Eq. 1",FaceColor="#d95f02"); hold on
h = histogram(C2,50,DisplayName="Eq. 2",FaceColor="#7570b3");
if runCase3 == true
histogram(C3,50,DisplayName="Eq. 2 + constant mixing",FaceColor="#e7298a");
end
if runCase4 == true
histogram(C4,50,DisplayName="Eq. 2 + time-varying mixing",FaceColor="#66a61e");
end
hold off

subplot(1,3,3)
plot(skLogn,kuLogn,DisplayName="Logn",Color="#1b9e77"); hold on
scatter(skewness(C1),kurtosis(C1),[],'MarkerEdgeColor','#000000','MarkerFaceColor','#d95f02');
scatter(skewness(C2),kurtosis(C2),[],'MarkerEdgeColor','#000000','MarkerFaceColor','#7570b3');
if runCase3 == true
    scatter(skewness(C3),kurtosis(C3),[],'MarkerEdgeColor','#000000','MarkerFaceColor','#e7298a');
end
if runCase4 == true
    scatter(skewness(C4),kurtosis(C4),[],'MarkerEdgeColor','#000000','MarkerFaceColor','#66a61e');
end
hold off
ylim([2 12]); xlim([0 3.5]);
legend('Log.','Eq. 1','Eq. 2','Eq. 2 + constant mixing','Eq. 2 + time-varying mixing',Location='southeast');

if runCase3 == true || runCase4 == true
    sgtitle("\mu_k = "+mu_k+", \sigma_k = "+sigma_k+", M = "+mix+"");
else
    sgtitle("\mu_k = "+mu_k+", \sigma_k = "+sigma_k+"");
end

%% Test basic solution of Case 1

mu = 0; t = 0:(1000/100000):1000;
sigma = 5e-4; eta = normrnd(0,1,length(t),1)';
X1 = exp(-mu*t);             % case one
X2 = exp(-(mu+sigma*eta).*t);  % case two

figure
plot(t,X1); hold on
plot(t,X2); hold off
disp(mean(X2));
legend("case one","case two");
title("\mu = " + mu + ", \sigma = " + sigma);

%% theoretical lognormal
% x = 1;
% 
% binEdgesOfConc = h.BinEdges;
% 
% X = binEdgesOfConc;
% 
% P = (1./(X.*sqrt(2*pi*sigma_k.^2)) ) .* exp(-(log(X)-mu_k).^2 / 2*sigma_k.^2);

% figure;plot(X,P);