function [] = plotKs(tr,ks,obs,sk,ku,obsLimA,obsLimB,EulLan,threshold,vuongRes,pV,limitOveride)
%plotKs
% INPUT: 
% OUTPUT: 

if nargin < 9
    threshold = 50;
end

if nargin < 8
    EulLan = true;
end

if nargin < 6
    obsLimA = 0.5;
    %obsLimB = length(tr);
    obsLimB = length(obs) + 0.5;
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

if nargin < 10
    vuongRes = 0;
end

if nargin == 12
    limits = limitOveride;
end

% Create Annotations for Vuong's Test Results
annot = strings(1,length(tr));
for i = 1:length(tr)
    if vuongRes(i) == 1
        annot(i) = "Normal";
    elseif vuongRes(i) == 2
        annot(i) = "Lognormal";
    elseif vuongRes(i) == 3
        annot(i) = "Weibull";
    elseif vuongRes(i) == 4
        annot(i) = "Gamma";
    elseif vuongRes(i) == 0
        annot(i) = "";
    end
end

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% Gamma family: generate theoretical skewness and kurtosis
kTh = linspace(0.2,5000,10000);
for i = 1:length(kTh)
    skGam(i) = 2/sqrt(kTh(i));
    kuGam(i) = 6/kTh(i) + 3;
end

% Weibull family: generate theoretical skewness and kurtosis
kWbl = linspace(0,5,10000);
for i = 1:length(kWbl)
    skWbl(i) = ( gamma(1 + 3/kWbl(i)) - 3*gamma(1 + 1/kWbl(i))*gamma(1 + 2/kWbl(i)) + 2*(gamma(1 + 1/kWbl(i)))^3 ) ./ ...
        ( gamma(1 + 2/kWbl(i)) -  (gamma(1 + 1/kWbl(i)))^2 )^(3/2);
    kuWbl(i) = ( gamma(1 + 4/kWbl(i)) - 4*gamma(1 + 1/kWbl(i))*gamma(1 + 3/kWbl(i)) + 6*( (gamma(1 + 1/kWbl(i)) )^2)*gamma(1 + 2/kWbl(i)) - 3*( (gamma(1 + 1/kWbl(i)))^4 ) ) ./ ...
       ( gamma(1 + 2/kWbl(i)) - ( gamma(1 + 1/kWbl(i)) )^2 )^2;
end

% Beta family: generate theoretical skewness and kurtosis
muB = linspace(0.5,5,100);
nuB = linspace(0.5,5,100);
alphaB = muB.*nuB; betaB = (1 - muB).*nuB;
for i = 1:length(betaB)
    skBet(i) = 2*(betaB(i) - alphaB(i)) * sqrt(alphaB(i) + betaB(i) + 1) ./ ...
        sqrt(betaB(i)*alphaB(i)) * (betaB(i) + alphaB(i) + 2);
    kuBet(i) = 3*(betaB(i) + alphaB(i) + 1) * ( 2*(betaB(i) + alphaB(i))^2 + alphaB(i)*betaB(i)*(betaB(i) + alphaB(i) - 6) ) ./ ...
        betaB(i)*alphaB(i)*(betaB(i) + alphaB(i) + 2)*(betaB(i) + alphaB(i) + 3);
end

% Loglogistic family: generate theoretical skewness and kurtosis
% Nothings shows up: Maybe I need to tune the shape parameter a bit (??)
% cLgi = linspace(2,4,10000);
% for i = 1:length(cLgi)
%     skLgi(i) = ( 2*(pi^2)*(csc(pi/cLgi(i)))^3 - 6*cLgi(i)*pi*csc(pi/cLgi(i))*csc(2*pi/cLgi(i)) + 3*(cLgi(i)^2)*csc(3*pi/cLgi(i)) ) ./ ...
%         sqrt(pi)*((-pi*(csc(pi/cLgi(i))^2)) + 2*cLgi(i)*csc(2*pi/cLgi(i)))^(3/2);
%     kuLgi(i) = ( -3*pi^2*(csc(pi/cLgi(i))^4) - 12*pi*cLgi(i)^2*csc(pi/cLgi(i))*csc(3*pi/cLgi(i)) + 4*cLgi(i)^3*csc(4*pi/cLgi(i)) + 6*cLgi(i)*pi^2*(csc(pi/cLgi(i))^3)*sec(pi/cLgi(i)) ) ./ ...
%         pi*( -pi*(csc(pi/cLgi(i))^2) + 2*cLgi(i)*csc(2*pi/cLgi(i)) )^2;
% end

% Only plot pV over 0.01
for i = 1:10
    for j = 1:length(tr)
        if pV(i,j) < 0.01
            pV(i,j) = nan;
        end
    end
end

% Only plot pV related to DOMINANT Vuong LLR
for k = 1:length(tr)
    if vuongRes(k) == 1
        disp('a...');
        pV([4 5 6 7 8 9 10],k) = nan;
    elseif vuongRes(k) == 2
        disp('b...');
        pV([2 3 4 7 8 9 10],k) = nan;
    elseif vuongRes(k) == 3
        disp('c..');
        pV([1 3 4 6 7 9 10],k) = nan;
    else
        disp('d...');
        pV([1 2 4 5 7 9 10],k) = nan;
    end
end

subplot(1,6,1)
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

subplot(1,6,[2 3])
plot(ks(1,:),tr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
hold on
plot(ks(2,:),tr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
plot(ks(3,:),tr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
plot(ks(4,:),tr,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
hold off
grid minor;
ylim(limits);
set(gca,'YDir','reverse');
legend(Location="best");
xlabel('p-value');
title('K-S p-values');

zzs = 0.25*ones(length(tr),1);
subplot(1,6,4)
% xline(0,HandleVisibility="off");
% hold on
% plot(pV(1,:),tr,DisplayName='Nor/Log',Color='#a6cee3',Marker='+',LineWidth=1);
% plot(pV(2,:),tr,DisplayName='Nor/Wbl',Color='#1f78b4',Marker='o',LineWidth=1);
% plot(pV(3,:),tr,DisplayName='Nor/Gam',Color='#b2df8a',Marker='*',LineWidth=1);
% plot(pV(5,:),tr,DisplayName='Log/Wbl',Color='#33a02c',Marker='.',LineWidth=1);
% plot(pV(6,:),tr,DisplayName='Log/Gam',Color='#fb9a99',Marker='x',LineWidth=1);
% plot(pV(8,:),tr,DisplayName='Wbl/Gam',Color='#e31a1c',Marker='square',LineWidth=1);
text(zzs,tr,annot,FontSize=8);
grid minor;
set(gca,"Xtick",[]);
ylim(limits); set(gca,'YDir','reverse');
% legend(Location="south");
title('Vuong LLR');

% subplot(1,6,4)
% plot(pV(1,:),tr,DisplayName='Nor/Log',Color='#a6cee3',Marker='+');
% hold on
% plot(pV(2,:),tr,DisplayName='Nor/Wbl',Color='#1f78b4',Marker='o');
% plot(pV(3,:),tr,DisplayName='Nor/Gam',Color='#b2df8a',Marker='*');
% plot(pV(5,:),tr,DisplayName='Log/Wbl',Color='#33a02c',Marker='.');
% plot(pV(6,:),tr,DisplayName='Log/Gam',Color='#fb9a99',Marker='x');
% plot(pV(8,:),tr,DisplayName='Wbl/Gam',Color='#e31a1c',Marker='square');
% hold off
% grid minor;
% ylim(limits); set(gca,'YDir','reverse');
% legend(Location="best");
% title('Vuong: p-values');

% subplot(1,5,4)
% yyaxis left
% plot(sk,tr,'DisplayName','Skewness'); hold on
% ylim(limits); set(gca,'YDir','reverse');
% % yticklabels({});
% yyaxis right
% plot(ku,tr,'DisplayName','Kurtosis');
% ylim(limits); set(gca,'YDir','reverse');
% % set(gca,'YTickLabel',{tr(5:5:length(tr))},'YColor','Black')
% xline(3,'.','Mesokurtic','HandleVisibility','off');
% xline(2.5,':','HandleVisibility','off');
% xline(3.5,':','HandleVisibility','off');
% xline(0,'.','Symmetric','HandleVisibility','off');
% xline(-0.5,':','HandleVisibility','off');
% xline(0.5,':','HandleVisibility','off');
% hold off
% grid minor;
% legend('Location','south');
% title('Moments');

subplot(1,6,[5 6])
numGroups = length(unique(tr));
clr = flipud(copper(numGroups));
gscatter(sk,ku,tr,clr);
hold on
% P = [0 2.5; 0 3.5; 0.5 3.5; 0.5 2.5] ;
% patch(P(:,1),P(:,2),[0.6 0.6 0.6],'EdgeColor','k','HandleVisibility','off');
plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#a6cee3',LineStyle='-',LineWidth=1);
plot(skGam,kuGam,'DisplayName','Gam.','Color','#1f78b4',LineStyle='--',LineWidth=1);
plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#b2df8a',LineStyle=':',LineWidth=1);
% plot(skBet,kuBet,'DisplayName','Beta','Color','#33a02c',LineStyle='-.',LineWidth=1);
% plot(skLgi,kuLgi,'DisplayName','Logl.','Color','#33a02c');
scatter(2,9,'DisplayName','Exp.',Marker='+',LineWidth=1);
scatter(0,9/5,'DisplayName','Uni.',Marker='o',LineWidth=1);
scatter(0,3,'DisplayName','Norm.',Marker='*',LineWidth=1);
scatter(0,21/5,'DisplayName','Logi.',Marker='.',LineWidth=1);
scatter(1.1395,5.4,'DisplayName','LEV',Marker='x',LineWidth=1);
hold off
grid minor;
ylim([1 10]); xlim([0 2.5]);
xlabel('Skewness'); ylabel('Kurtosis');
lgd = legend('Location','best','FontSize',6);
title(lgd,'P [dbar]');
title('Skewness vs. Kurtosis');

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