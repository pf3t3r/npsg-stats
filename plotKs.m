function [] = plotKs(tr,ks,obs,sk,ku,obsLimA,obsLimB,EulLan,threshold,limitOveride)
%plotKs
% INPUT: 
% OUTPUT: 

if nargin < 9
    threshold = 100;
end

if nargin < 8
    EulLan = true;
end

if nargin < 6
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

if nargin == 10
    limits = limitOveride;
end

subplot(1,3,1)
barh(obs,'FaceColor','#a6cee3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
set(gca,'XDir','reverse');
ylim([obsLimA obsLimB]);
ylabel('Pressure [dbar]');
set(gca,"YTick",1:1:length(ytix),"YTickLabel",ytix);
title('No. of Observations');

subplot(1,3,2)
plot(ks(1,:),tr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
hold on
plot(ks(2,:),tr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
plot(ks(3,:),tr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
plot(ks(4,:),tr,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
hold off
grid minor;
ylim(limits);
set(gca,'YDir','reverse');
legend('Location','best');
xlabel('p-value');
% ylabel('Pressure [db]');
title('KS Test');

subplot(1,3,3)
yyaxis left
plot(sk,tr,'DisplayName','Skewness'); hold on
ylim(limits); set(gca,'YDir','reverse');
% yticklabels({});
yyaxis right
plot(ku,tr,'DisplayName','Kurtosis');
ylim(limits); set(gca,'YDir','reverse');
% set(gca,'YTickLabel',{tr(5:5:length(tr))},'YColor','Black')
xline(3,'.','Mesokurtic','HandleVisibility','off');
xline(2.5,':','HandleVisibility','off');
xline(3.5,':','HandleVisibility','off');
xline(0,'.','Symmetric','HandleVisibility','off');
xline(-0.5,':','HandleVisibility','off');
xline(0.5,':','HandleVisibility','off');
hold off
grid minor;
legend('Location','south');
title('Moments');

% old subplot 3
% subplot(1,3,3)
% plot(sk,tr,'Color','#1f78b4');
% ylim(limits);
% tAx = gca;
% xlabel(tAx,'Skewness','Color','#1f78b4');
% % tAx2 = axes('Position', get(tAx, 'Position')); % Create a new axes in the same position as the first one, overlaid on top
% tAx2 = axes('Position', [0.6916    0.1060    0.2134    0.7825]); % Create a new axes in the same position as the first one, overlaid on top
% plot(ku,tr,'Color','#33a02c'); 
% set(tAx2, 'ylim', get(tAx, 'ylim'), 'color', 'none'); % Set y limits same as original axes, and make background transparent
% set(tAx2,'YTickLabel',[]);
% set(tAx,'YDir','reverse');
% set(tAx2,'YDir','reverse');
% xlabel(tAx2,'Kurtosis','Color','#33a02c');
% % tAx2.set('XTickLabel','Color','#33a02c');
% set(tAx2,'XAxisLocation','top');
% set(tAx,'XColor','#1f78b4');
% set(tAx2,'XColor','#33a02c');

% data1m=get(tAx,'ylim');
% data2m=get(tAx2,'ylim');
% nt = 10; %tick number
% data1_tick=linspace(data1m(1),data1m(2),nt);
% data2_tick=linspace(data2m(1),data2m(2),nt);
% set(tAx,'ytick',round(data1_tick*100)/100);
% set(tAx2,'ytick',round(data2_tick*100)/100);

% subplot(1,4,3)
% plot(sk,tr,'Color','#a6cee3','LineWidth',2);
% set(gca,'YDir','reverse');
% title('Skewness');
% ylim(limits);
% 
% subplot(1,4,4)
% plot(ku,tr,'Color','#1f78b4','LineWidth',2);
% set(gca,'YDir','reverse');
% title('Kurtosis');
% ylim(limits);

end