% File to open and examine HOTS chloropigment data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 29 15]);
set(0,'defaultAxesFontSize',12);

addpath("baroneRoutines\");

%% Try Barone's Functions: bbvuong

% load("datafiles\chloro.mat","chloro200",'pgrid200');
% 
% [testR,testp2] = bbvuong(chloro200(1:100,1));
% normal dist fits best here....

% [testR_t,testp2_t] = bbvuong(chloro2D(50,:));

% there is an MLE function in MATLAB!
% [phat,pci] = mle(chloro2D(1:100,1));

%% Try: hotFTPextract

% cruise = 300;
% filename = 'datafiles\ctd_iso' + string(cruise) + '.mat';
% 
% if isfile(filename)
%     disp('CTD and ISO already created for that cruise');
% else
%     disp('Extracting...');
%     [ctd,iso] = hotFTPextract(cruise);
%     save(filename,"ctd","iso");
% end

%% Can I save multiple cruises easily?

% process broke after cruise 240 so extract 241-329 now... merge later
% filename = 'datafiles/ctd_iso_master2';
% 
cruises = 329;
missingCruises = [21,207];
cruisesRecorded = linspace(1,cruises,cruises);
cruisesRecorded(missingCruises) = [];
clear cruises missingCruises;
% 
% cruisesRemaining = cruisesRecorded(239:end);

%% don't run again

% dontRun = 0;
% 
% if dontRun == 0
%     for i = cruisesRemaining
%         [ctd(i),iso(i)] = hotFTPextract(i);
%         save(filename,"ctd","iso");
%         disp(i);
%     end
% else
%     disp('already calculated..');
% end

% No need to run this again!

% save('datafiles/ctd_iso_ALL','ctd','iso');

%% So let's extract the fluorescence and group them by casts
% also let's extract the pcm

load('datafiles\ctd_iso_ALL.mat');

f = [];
t = [];
pcm = [];
for i = cruisesRecorded
    tmpF = ctd(i).f;
    tmpT = ctd(i).date;
    tmpPCM = ctd(i).pcm;
    f = [f tmpF];
    t = [t tmpT];
    pcm = [pcm tmpPCM];
end
clear tmpF tmpT tmpPCM i;
t = datetime(t,'ConvertFrom','datenum');

f(f<0) = NaN;
depthMeasurements = 129;
eulerianDepth = linspace(0,2*depthMeasurements,depthMeasurements);
lagrangianDepth = linspace(-128,128,depthMeasurements);

%%
% figure
% plot(f_e(:,3109:3112),eulerianDepth,'b-','DisplayName','June 26 2008');
% hold on
% plot(f_e(:,3113:3121),eulerianDepth,'r--','DisplayName','July 26 2008');
% hold off
% legend();
% set(gca,"YDir","reverse");

%% remove casts for which the pressure of the CM is NAN
% (where the CM could not be found)

f_copy = f;

p_copy = pcm;
p_copy(isnan(pcm)) = [];

f_copy(:,isnan(pcm)) = [];
f_e = f_copy(1:129,:);

t_copy = t;
t_copy(isnan(pcm)) = [];

%% flag (visually for now) casts with abnormally large values
clc;
for i = 1:4380
%     if mean(f_e(:,i)) > 0.4
%         disp(i);
%     end
    for j = 1:129
        if f_e(j,i) > 2
            disp(i);
        end
    end
end

flagged = f_e(:,[2434 2599 2626 2898 2918 2926]);
save('datafiles\flagged.mat',"flagged");
%% Remove flagged (visually confirmed) casts

dodgyCasts = [18 248 253 666 667 668 669 ...
    1213 1214 1215 1216 1217 1218 1219 1220 1221 1222 1223 1224 ...
    2124 2371 2434 2599 2610 2626 2898 2918 2926 2979 ...
    3110 3133 3560 3561];
f_er = f_e;
f_er(:,dodgyCasts) = [];
t_copy(dodgyCasts) = [];
save('datafiles\f_all.mat',"f_er","t_copy");

%%
figure
plot(f_er,eulerianDepth,'Color',[0.6 0.6 0.6]);
set(gca,"YDir","reverse");

%% Convert Eulerian f into Lagrangian coordinates

% move lag stuff to after removal of flagged data
f_lag = NaN(129,length(f_er(1,:)));
offset = floor(0.5*(129 - p_copy));

%%
for i = 1:length(f_er(1,:))
    f_lag(:,i) = circshift(f_er(:,i),offset(i));
    if offset(i) > -1 && offset(i) < 40
        disp(i);
        f_lag(1:offset(i),i) = NaN;
    elseif offset(i) == -1
        f_lag(end,i) = NaN;
    elseif offset(i) < -1 && offset(i) > -40
        disp(i);
        f_lag((end+offset(i)):end,i) = NaN;
    elseif abs(offset(i)) > 40
        f_lag(:,i) = NaN;
    end
end

save('datafiles\f_all.mat',"f_lag",'-append');

save('datafiles\f_all.mat',"t_copy",'-append');

%%
figure
plot(f_lag,lagrangianDepth);

%% check crn and cast of suspicious casts
% load('C:\Users\pfarrell\OneDrive - NIOZ\Documenten\GitHub\PhD-chloro\datafiles\ctd_iso_master.mat');
% 
% test = datetime(ctd(10).date(1),'ConvertFrom','datenum');

%% Let's try KS Test on f

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = f_er(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksE(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = f_lag(i,:);
    tmp(isnan(tmp)) = [];
    [~,ksL(:,i),~] = statsplot2(tmp,'noplot');
end
clear tmp;

%% Plot the above
ax = figure;
ax.Position = [3 3 20 15];

% Eulerian
subplot(1,2,1)
plot(ksE(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksE(2,:),eulerianDepth,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksE(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE(4,:),eulerianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend();
ylim([0 250]);
set(gca,'YDir','reverse');
xlabel('p-value');
ylabel('Depth [m]');
title('Eulerian');

% Lagrangian
subplot(1,2,2)
plot(ksL(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL(2,:),lagrangianDepth,'+--','Color',[0 0 0],'LineStyle','--','DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL(4,:),lagrangianDepth,'r.--','DisplayName','Gamma','MarkerSize',4);
hold off
legend();
ylim([-125 125]);
set(gca,'YDir','reverse');
xlabel('p-value');
ylabel('Depth [m]');
title('Lagrangian');

sgtitle('Kolmogorov-Smirnov Test (1989-2021)');
exportgraphics(ax,'figures/ks_allCast_89-21.png');

%% histfit practice

% for i = 1:129
%     text = 'lvl' + string(i);
%     savefileName = 'figures/lvlBylvl/hist_' + text + '.png';
% 
%     tmp = f_er(i,:);
%     tmp(isnan(tmp)) = [];
% 
%     text = figure;
%     [~,ksE(:,i),~] = statsplot2_pf(tmp,100,i);
% 
%     exportgraphics(text,savefileName);
% end

%%

% test = makedist("Normal",mean(tmp),std(tmp));
pd = makedist('Normal','mu',mean(tmp),'sigma',std(tmp));
pd2 = makedist('Lognormal','mu',mean(log(tmp)),'sigma',std(log(tmp)));
% pd3 = makedist('Weibull',);
% x = (mean(tmp)-5*std(tmp)):(mean(tmp)+5*std(tmp));
% y = pdf(pd, x);

figure;
% histogram(tmp);
% hold on
yyaxis left;
plot(pd);
hold on
plot(pd2);
yyaxis right;
histogram(tmp);
hold off
% hold off
legend();


%% isopyncnal view of DCM
% plot sig vs. f

% load(filename);
% 
% % mean fluorescence for that cruise
% f_mean = mean(ctd.f,2,'omitnan');
% sig_mean = mean(ctd.sig,2,"omitnan");
% 
% figureName = 'figures/isopycnal_DCM_' + string(cruise) + '.png';
% 
% ax1 = figure;
% scatter(ctd.f,-ctd.sig,'Marker','.','MarkerEdgeColor',[0.6 0.6 0.6],'HandleVisibility','off');
% hold on
% plot(f_mean,-sig_mean,'Color','red','LineWidth',1.5,'DisplayName','Mean Fluorescence');
% hold off
% legend();
% xlabel('Chloropigment [ug L^{-1}]');
% ylabel('Potential Density Anomaly \sigma_0 [kg m^{-3}]');
% title(sprintf('Isopycnal View of DCM: \\sigma_0 vs fluorescence (Cruise %s)',string(cruise)));
% exportgraphics(ax1,figureName);

%% Try: RunMedian

% testDataRun = RunMedian(chloro200(:,1),11);

%% plot the above vs original
% figure
% plot(chloro200(:,1),-pgrid200(:,1));
% hold on
% plot(testDataRun,-pgrid200(:,1));
% hold off

%% Try: statsplot2

% [MLEp,KSp,nll] = statsplot2(chloro200(1:100,1));
% clear MLEp KSp nll

%% Try kstest

% h_norm = kstest(chloro200(:,1));
% if testh = 1, then it means that the distribution is normal (?)
