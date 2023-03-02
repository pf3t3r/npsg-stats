% File to open and examine HOTS chloropigment data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',16);

addpath("baroneRoutines\");

%% Try Barone's Functions: bbvuong

load("datafiles\chloro.mat","chloro200",'pgrid200');

[testR,testp2] = bbvuong(chloro200(1:100,1));
% normal dist fits best here....

% [testR_t,testp2_t] = bbvuong(chloro2D(50,:));

% there is an MLE function in MATLAB!
% [phat,pci] = mle(chloro2D(1:100,1));

%% Try: hotFTPextract

cruise = 300;
filename = 'datafiles\ctd_iso' + string(cruise) + '.mat';

if isfile(filename)
    disp('CTD and ISO already created for that cruise');
else
    disp('Extracting...');
    [ctd,iso] = hotFTPextract(cruise);
    save(filename,"ctd","iso");
end

%% Can I save multiple cruises easily?
filename = 'datafiles/ctd_iso_master';

cruises = 329;
missingCruises = [21,207];
cruisesRecorded = linspace(1,cruises,cruises);
cruisesRecorded(missingCruises) = [];

for i = cruisesRecorded
    [ctd(i),iso(i)] = hotFTPextract(i);
    save(filename,"ctd","iso");
    disp(i);
end

% No need to run this again!

%% So let's extract the fluorescence and group them by casts
% also let's extract the pcm

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

%% Convert Eulerian f into Lagrangian coordinates

% First remove columns where CM = NaN
f_copy = f;
p_copy = pcm;
f_copy(:,isnan(pcm)) = [];
p_copy(isnan(pcm)) = [];

f_e = f_copy(1:129,:);
f_lag = NaN(129,length(p_copy));
offset = floor(0.5*(129 - p_copy));

%%
for i = 1:length(p_copy)
    f_lag(:,i) = circshift(f_e(:,i),offset(i));
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

%%
figure
plot(f_e(:,1224),eulerianDepth);
set(gca,"YDir","reverse");

%% test for dodgy fluorescence
clc;
for i = 1:4380
%     if mean(f_e(:,i)) > 0.4
%         disp(i);
%     end
    for j = 2:129
        if f_e(j,i) > 2
            disp(i);
        end
    end
end

%% Eulerian Fluorescence with dodgy casts removed

dodgyCasts = [18 248 253 666 667 668 669 ...
    1213 1214 1215 1216 1217 1218 1219 1220 1221 1222 1223 1224 ...
    2124 2371 2610 2979 ...
    3110 3133 3560 3561];
f_er = f_e;
f_er(:,dodgyCasts) = [];

%%
figure
plot(f_er,eulerianDepth,'Color',[0.6 0.6 0.6]);
set(gca,"YDir","reverse");
%%
figure
plot(f_lag(:,1:4000),lagrangianDepth);

%% Let's try KS Test on f

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = f(i,:);
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
% Eulerian
% subplot(1,2,1)
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
title('Eulerian KS-Test (89-21)');
ax.Position = [3 3 20 15];

%% isopyncnal view of DCM
% plot sig vs. f

load(filename);

% mean fluorescence for that cruise
f_mean = mean(ctd.f,2,'omitnan');
sig_mean = mean(ctd.sig,2,"omitnan");

figureName = 'figures/isopycnal_DCM_' + string(cruise) + '.png';

ax1 = figure;
scatter(ctd.f,-ctd.sig,'Marker','.','MarkerEdgeColor',[0.6 0.6 0.6],'HandleVisibility','off');
hold on
plot(f_mean,-sig_mean,'Color','red','LineWidth',1.5,'DisplayName','Mean Fluorescence');
hold off
legend();
xlabel('Chloropigment [ug L^{-1}]');
ylabel('Potential Density Anomaly \sigma_0 [kg m^{-3}]');
title(sprintf('Isopycnal View of DCM: \\sigma_0 vs fluorescence (Cruise %s)',string(cruise)));
exportgraphics(ax1,figureName);

%% Try: RunMedian

testDataRun = RunMedian(chloro200(:,1),11);

%% plot the above vs original
figure
plot(chloro200(:,1),-pgrid200(:,1));
hold on
plot(testDataRun,-pgrid200(:,1));
hold off

%% Try: statsplot2

% [MLEp,KSp,nll] = statsplot2(chloro200(1:100,1));
% clear MLEp KSp nll

%% Try kstest

h_norm = kstest(chloro200(:,1));
% if testh = 1, then it means that the distribution is normal (?)
