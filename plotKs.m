function [] = plotKs(tr,ks,obs,sk,obsLimA,obsLimB,EulLan,threshold)
%plotKs
% INPUT: 
% OUTPUT: 

if nargin < 8
    threshold = 100;
end

if nargin < 7
    EulLan = true;
end

if nargin < 5
    obsLimA = 1;
    obsLimB = length(tr);
end

if EulLan
    limits = [0 200];
    if obsLimB-obsLimA > 20
        ytix = 5:5:205;
    else
        ytix = 5:10:205;
    end
else
    limits = [-150 100];
    ytix = tr;
end

subplot(1,3,1)
barh(obs,'FaceColor','#a6cee3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([obsLimA obsLimB]);
ylabel('Pressure [db]');
set(gca,"YTick",1:1:length(ytix),"YTickLabel",ytix);
title('No. of Observations');

subplot(1,3,2)
plot(ks(1,:),tr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
hold on
plot(ks(2,:),tr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
plot(ks(3,:),tr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
plot(ks(4,:),tr,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
% plot(ks(1,:),tr,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
% hold on
% plot(ks(2,:),tr,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
% plot(ks(3,:),tr,'xr-','DisplayName','Weibull','MarkerSize',4);
% plot(ks(4,:),tr,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
grid minor;
ylim(limits);
% xlim([0 0.9]);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
% ylabel('Pressure [db]');
title('KS Test');

subplot(1,3,3)
plot(sk,tr,'Color','#a6cee3','LineWidth',2);
set(gca,'YDir','reverse');
title('Skewness');
ylim(limits);

end