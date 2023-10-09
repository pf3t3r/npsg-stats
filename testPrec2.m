clc; close all; clear;

%% What is real mean and SD of chl in L1?

% tmp = importdata('data\L1\hplcChla_88-21_150.txt').data(:,5);
% 
% testMean = mean(tmp);
% testSd = std(tmp);

%% Initialise MEAN and STD

mu = [0.5 1 1.5 2 2.5 3];
sd = 0.4;

% mu = round(testMean); sd = round(testSd);
% Lognormal
% rl20 = nan(20,5); rl40 = nan(40,5); rl60 = nan(60,5); rl80 = nan(80,5); rl100 = nan(100,5); 
% rl20(:,1) = lognrnd(mu,sd,20,1); rl40(:,1) = lognrnd(mu,sd,40,1);
% rl60(:,1) = lognrnd(mu,sd,60,1); rl80(:,1) = lognrnd(mu,sd,80,1);
% rl100(:,1) = lognrnd(mu,sd,100,1);
% 
% % Normal
% rn20 = nan(20,5); rn40 = nan(40,5); rn60 = nan(60,5); rn80 = nan(80,5); rn100 = nan(100,5);
% 
% % Use rand or...
% % rn20(:,1) = rand(20,1); rn40(:,1) = rand(40,1); rn60(:,1) = rand(60,1);
% % rn80(:,1) = rand(80,1); rn100(:,1) = rand(100,1);
% 
% % ...normrnd
% rn20(:,1) = abs(normrnd(mu,sd,20,1)); rn40(:,1) = abs(normrnd(mu,sd,40,1));
% rn60(:,1) = abs(normrnd(mu,sd,60,1)); rn80(:,1) = abs(normrnd(mu,sd,80,1));
% rn100(:,1) = abs(normrnd(mu,sd,100,1));

%%

% figure;
% subplot(2,1,1)
% histogram(rl40);
% subplot(2,1,2)
% histogram(log(rl40));
%% Rounding

% % Lognormal
% [rl20(:,2),rl20(:,3),rl20(:,4),rl20(:,5)] = quickRound(rl20);
% [rl40(:,2),rl40(:,3),rl40(:,4),rl40(:,5)] = quickRound(rl40);
% [rl60(:,2),rl60(:,3),rl60(:,4),rl60(:,5)] = quickRound(rl60);
% [rl80(:,2),rl80(:,3),rl80(:,4),rl80(:,5)] = quickRound(rl80);
% [rl100(:,2),rl100(:,3),rl100(:,4),rl100(:,5)] = quickRound(rl100);
% 
% % Normal
% [rn20(:,2),rn20(:,3),rn20(:,4),rn20(:,5)] = quickRound(rn20);
% [rn40(:,2),rn40(:,3),rn40(:,4),rn40(:,5)] = quickRound(rn40);
% [rn60(:,2),rn60(:,3),rn60(:,4),rn60(:,5)] = quickRound(rn60);
% [rn80(:,2),rn80(:,3),rn80(:,4),rn80(:,5)] = quickRound(rn80);
% [rn100(:,2),rn100(:,3),rn100(:,4),rn100(:,5)] = quickRound(rn100);
% 
% rn20(rn20<0) = 0.1; rn40(rn40<0) = 0.1; rn60(rn60<0) = 0.1; rn80(rn80<0) = 0.1;
%% test

% % log on log
% [ksL20,lil_LL20] = quickMleAndKs(rl20);
% ksL40 = quickMleAndKs(rl40);
% ksL60 = quickMleAndKs(rl60);
% ksL80 = quickMleAndKs(rl80);
% ksL100 = quickMleAndKs(rl100);
% 
% % norm on norm
% ksN20 = quickMleAndKs(rn20,'norm');
% ksN40 = quickMleAndKs(rn40,'norm');
% ksN60 = quickMleAndKs(rn60,'norm');
% ksN80 = quickMleAndKs(rn80,'norm');
% ksN100 = quickMleAndKs(rn100,'norm');
% 
% % log on norm
% ks_ln_20 = quickMleAndKs(rn20);
% ks_ln_40 = quickMleAndKs(rn40);
% ks_ln_60 = quickMleAndKs(rn60);
% ks_ln_80 = quickMleAndKs(rn80);
% ks_ln_100 = quickMleAndKs(rn100);
% 
% % norm on long
% ks_nl_20 = quickMleAndKs(rl20,'norm');
% ks_nl_40 = quickMleAndKs(rl40,'norm');
% ks_nl_60 = quickMleAndKs(rl60,'norm');
% ks_nl_80 = quickMleAndKs(rl80,'norm');
% ks_nl_100 = quickMleAndKs(rl100,'norm');

%% MLE
% for i = 1:5
%     mle20(i,:) = mle(rl20(:,i),'distribution','Lognormal');
%     mle40(i,:) = mle(rl40(:,i),'distribution','Lognormal');
%     mle60(i,:) = mle(rl60(:,i),'distribution','Lognormal');
%     mle80(i,:) = mle(rl80(:,i),'distribution','Lognormal');
%     mle100(i,:) = mle(rl100(:,i),'distribution','Lognormal');
% 
%     xCdf20(i,:) = linspace(min(rl20(:,i)) - 2*std(rl20(:,i)), max(rl20(:,i)) + 2*std(rl20(:,i)),2000);
%     xCdf40(i,:) = linspace(min(rl40(:,i)) - 2*std(rl40(:,i)), max(rl40(:,i)) + 2*std(rl40(:,i)),2000);
%     xCdf60(i,:) = linspace(min(rl60(:,i)) - 2*std(rl60(:,i)), max(rl60(:,i)) + 2*std(rl60(:,i)),2000);
%     xCdf80(i,:) = linspace(min(rl80(:,i)) - 2*std(rl80(:,i)), max(rl80(:,i)) + 2*std(rl80(:,i)),2000);
%     xCdf100(i,:) = linspace(min(rl100(:,i)) - 2*std(rl100(:,i)), max(rl100(:,i)) + 2*std(rl100(:,i)),2000);
% 
%     yCdfLogn20(i,:) = cdf('logn', xCdf20(i,:), mle20(i,1), mle20(i,2));
%     yCdfLogn40(i,:) = cdf('logn', xCdf40(i,:), mle40(i,1), mle40(i,2));
%     yCdfLogn60(i,:) = cdf('logn', xCdf60(i,:), mle60(i,1), mle60(i,2));
%     yCdfLogn80(i,:) = cdf('logn', xCdf80(i,:), mle80(i,1), mle80(i,2));
%     yCdfLogn100(i,:) = cdf('logn', xCdf100(i,:), mle100(i,1), mle100(i,2));
% 
% end
% KS
% for i = 1:5
%     [Hl_20(i),KS_20(i)] = kstest(rl20(:,i),[xCdf20(i,:)' yCdfLogn20(i,:)']);
%     [Hl_40(i),KS_40(i)] = kstest(rl40(:,i),[xCdf40(i,:)' yCdfLogn40(i,:)']);
%     [Hl_60(i),KS_60(i)] = kstest(rl60(:,i),[xCdf60(i,:)' yCdfLogn60(i,:)']);
%     [Hl_80(i),KS_80(i)] = kstest(rl80(:,i),[xCdf80(i,:)' yCdfLogn80(i,:)']);
%     [Hl_100(i),KS_100(i)] = kstest(rl100(:,i),[xCdf100(i,:)' yCdfLogn100(i,:)']);
% end

%%
% noSig = [5 4 3 2 1];
% 
% ax1 = figure;
% subplot(2,2,1)
% plot(noSig,ksL20,'DisplayName','20');
% hold on
% plot(noSig,ksL40,'DisplayName','40');
% plot(noSig,ksL60,'DisplayName','60');
% plot(noSig,ksL80,'DisplayName','80');
% plot(noSig,ksL100,'DisplayName','100');
% hold off
% grid on;
% xlabel('No. of significant digits'); xlim([1 4]);
% ylabel('p-value'); ylim([0 1]);
% title("Logn test of logn data ($\mu$ = " + mu + ", $\sigma$ = " + sd + ")",'Interpreter','latex');
% 
% subplot(2,2,2)
% plot(noSig,ksN20,'DisplayName','20');
% hold on
% plot(noSig,ksN40,'DisplayName','40');
% plot(noSig,ksN60,'DisplayName','60');
% plot(noSig,ksN80,'DisplayName','80');
% plot(noSig,ksN100,'DisplayName','100');
% hold off
% grid on;
% xlabel('No. of significant digits'); xlim([1 4]);
% ylabel('p-value'); ylim([0 1]);
% title("Norm test of norm data ($\mu$ = " + mu + ", $\sigma$ = " + sd + ")",'Interpreter','latex');
% 
% subplot(2,2,3)
% plot(noSig,ks_nl_20,'DisplayName','20');
% hold on
% plot(noSig,ks_nl_40,'DisplayName','40');
% plot(noSig,ks_nl_60,'DisplayName','60');
% plot(noSig,ks_nl_80,'DisplayName','80');
% plot(noSig,ks_nl_100,'DisplayName','100');
% hold off
% grid on;
% xlabel('No. of significant digits'); xlim([1 4]);
% ylabel('p-value'); ylim([0 1]);
% title("Norm test of logn data ($\mu$ = " + mu + ", $\sigma$ = " + sd + ")",'Interpreter','latex');
% lgd = legend('Location','best');
% lgd.Title.String = 'No. of samples';
% 
% subplot(2,2,4)
% plot(noSig,ks_ln_20,'DisplayName','20');
% hold on
% plot(noSig,ks_ln_40,'DisplayName','40');
% plot(noSig,ks_ln_60,'DisplayName','60');
% plot(noSig,ks_ln_80,'DisplayName','80');
% plot(noSig,ks_ln_100,'DisplayName','100');
% hold off
% grid on;
% xlabel('No. of significant digits'); xlim([1 4]);
% ylabel('p-value'); ylim([0 1]);
% title("Logn test of norm data ($\mu$ = " + mu + ", $\sigma$ = " + sd + ")",'Interpreter','latex');
% 
% sgtitle('KS Precision Sensitivity');
% exportgraphics(ax1,'figures/prec/_mu6_s1.png');
% 
% clear i;

%% Restart:

% run 50x, then take mean

tmp20 = []; tmp40 = []; tmp60 = []; tmp80 = []; tmp100 = []; tmp1000 = [];

for j = 1:length(mu)
for i = 1:50
    t20 = lognrnd(mu(j),sd,20,1);
    tmp20 = [tmp20 t20];
    t40 = lognrnd(mu(j),sd,40,1);
    tmp40 = [tmp40 t40];
    t60 = lognrnd(mu(j),sd,60,1);
    tmp60 = [tmp60 t60];
    t80 = lognrnd(mu(j),sd,80,1);
    tmp80 = [tmp80 t80];
    t100 = lognrnd(mu(j),sd,100,1);
    tmp100 = [tmp100 t100];
    t1000 = lognrnd(mu(j),sd,1000,1);
    tmp1000 = [tmp1000 t1000];
end
r2l20(:,j) = mean(tmp20,2); r2l40(:,j) = mean(tmp40,2);
r2l60(:,j) = mean(tmp60,2); r2l80(:,j) = mean(tmp80,2);
r2l100(:,j) = mean(tmp100,2); r2l1000(:,j) = mean(tmp1000,2);
end 

figure;
% plot(tmp20(:,1:50),'g:');
plot(r2l20);
hold on
plot(r2l40);
plot(r2l60);
plot(r2l80);
plot(r2l100);
plot(r2l1000);
hold off
legend();

%%

% log on log
ksL20 = quickMleAndKs(r2l20,6);
ksL40 = quickMleAndKs(r2l40,6);
ksL60 = quickMleAndKs(r2l60,6);
ksL80 = quickMleAndKs(r2l80,6);
ksL100 = quickMleAndKs(r2l100,6);
ksL1000 = quickMleAndKs(r2l1000,6);

% save("output\Prec\mu6.mat","ksL20","ksL40","ksL60","ksL80","ksL100","ksL1000");

%%

% mu = [2 3 4 5 6];
% % 20 values
% tmp1 = load("output\Prec\mu2.mat").ksL20;
% tmp2 = load("output\Prec\mu3.mat").ksL20;
% tmp3 = load("output\Prec\mu4.mat").ksL20;
% tmp4 = load("output\Prec\mu5.mat").ksL20;
% tmp5 = load("output\Prec\mu6.mat").ksL20;
% ttmp20 = [tmp1 tmp2 tmp3 tmp4 tmp5];
% 
% % 40 values
% tmp1 = load("output\Prec\mu2.mat").ksL40;
% tmp2 = load("output\Prec\mu3.mat").ksL40;
% tmp3 = load("output\Prec\mu4.mat").ksL40;
% tmp4 = load("output\Prec\mu5.mat").ksL40;
% tmp5 = load("output\Prec\mu6.mat").ksL40;
% ttmp40 = [tmp1 tmp2 tmp3 tmp4 tmp5];
% 
% % 60 values
% tmp1 = load("output\Prec\mu2.mat").ksL60;
% tmp2 = load("output\Prec\mu3.mat").ksL60;
% tmp3 = load("output\Prec\mu4.mat").ksL60;
% tmp4 = load("output\Prec\mu5.mat").ksL60;
% tmp5 = load("output\Prec\mu6.mat").ksL60;
% ttmp60 = [tmp1 tmp2 tmp3 tmp4 tmp5];
% 
% % 80 values
% tmp1 = load("output\Prec\mu2.mat").ksL80;
% tmp2 = load("output\Prec\mu3.mat").ksL80;
% tmp3 = load("output\Prec\mu4.mat").ksL80;
% tmp4 = load("output\Prec\mu5.mat").ksL80;
% tmp5 = load("output\Prec\mu6.mat").ksL80;
% ttmp80 = [tmp1 tmp2 tmp3 tmp4 tmp5];
% 
% % 100 values
% tmp1 = load("output\Prec\mu2.mat").ksL100;
% tmp2 = load("output\Prec\mu3.mat").ksL100;
% tmp3 = load("output\Prec\mu4.mat").ksL100;
% tmp4 = load("output\Prec\mu5.mat").ksL100;
% tmp5 = load("output\Prec\mu6.mat").ksL100;
% ttmp100 = [tmp1 tmp2 tmp3 tmp4 tmp5];
% 
% % 1000 values
% tmp1 = load("output\Prec\mu2.mat").ksL1000;
% tmp2 = load("output\Prec\mu3.mat").ksL1000;
% tmp3 = load("output\Prec\mu4.mat").ksL1000;
% tmp4 = load("output\Prec\mu5.mat").ksL1000;
% tmp5 = load("output\Prec\mu6.mat").ksL1000;
% ttmp1000 = [tmp1 tmp2 tmp3 tmp4 tmp5];

ax = figure;
plot(mu,ksL20,'DisplayName','n = 20');
hold on
plot(mu,ksL40,'DisplayName','n = 40');
plot(mu,ksL60,'DisplayName','n = 60');
plot(mu,ksL80,'DisplayName','n = 80');
plot(mu,ksL100,'DisplayName','n = 100');
plot(mu,ksL1000,'DisplayName','n = 1000');
hold off
ylabel('p-value');
xlabel('mean');
legend();
title('K-S Lognormal Test of Random Lognormal Data','Sensitivity Analysis (\sigma = 1)');
exportgraphics(ax,'figures/prec/run5.png'); clear ax;