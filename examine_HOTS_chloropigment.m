% File to open and examine HOTS chloropigment data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',16);

%% Open data file and extract variables
data = importdata('data/hots-chloropigment.txt').data;

crn = data(:,1);
day = data(:,2);
pressure = data(:,3);
chloro = data(:,4);

chloro(chloro==-9) = NaN;

%% Examine First Days of Data

ax1 = figure;
plot(chloro(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
hold on
plot(chloro(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
plot(chloro(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
hold off
legend();
ylabel('pressure [db]');
xlabel('chloropigment (fluorescence) [ug/L]');
title('CTD: chloropigment vs depth (HOTS 88/89)');

exportgraphics(ax1,'figures/ctd-day-1.png');

%% Fluorescence Depth- and Time-Series (1988-2022): Eulerian View

chloro2D = reshape(chloro,101,[]);
days = day + datetime(1988,09,30);
days = reshape(days,101,[]);
pres = reshape(pressure,101,[]);

[t_grid,p_grid] = meshgrid(datenum(days(1,:)),pres(:,1));

ax2 = figure;
contourf(t_grid,p_grid,chloro2D,0:0.14:1.4,'LineColor','auto');
set(gca,'Ydir','reverse')
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
xlabel('Time');
ylabel('Depth [db]');
title('DCM Time Series: 1988 - 2021 (Eulerian Perspective)');

exportgraphics(ax2,'figures/fluorescence-1988-2021_eulerianView.png');

%% Fluorescence Depth- and Time-Series (1988-2022): Lagrangian View

% Find the DCM
for i=1:329
    [val(i),idx(i)] = max(chloro2D(:,i));
end

% Put pressure in terms of DCM
for i = 1:329
    p_lang(:,i) = pres(:,i) - pres(idx(i),i);
end

% Put fluorescence in terms of DCM
% Shift the original fluorescence data such that the DCM is centred
chloro_lang = zeros(101,329);
midpt = 51;
offset = midpt - idx;

for i = 1:329
    chloro_lang(:,i) = circshift(chloro2D(:,i),offset(i));
    if offset(i) > -1 && offset(i) < 40
        disp(i);
        chloro_lang(1:offset(i),i) = NaN;
    elseif offset(i) == -1
        chloro_lang(end,i) = NaN;
    elseif offset(i) < -1 && offset(i) > -40
        disp(i);
        chloro_lang((end+offset(i)):end,i) = NaN;
    elseif abs(offset(i)) > 40
        chloro_lang(:,i) = NaN;
    end
end

% Create meshgrid for time and pressure in Lagrangian view
[t_lang_grid,p_lang_grid] = meshgrid(datenum(days(1,:)),pres(:,1)-100);

% Make a filled contour plot of the DCM in the Lagrangian perspective
ax3 = figure;
contourf(t_lang_grid,p_lang_grid,chloro_lang,0:0.14:1.4,'LineColor','auto');
set(gca,'Ydir','reverse')
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
xlabel('Time');
ylabel('Depth [db]');
title('DCM Time Series: 1988 - 2021 (Lagrangian Perspective)');

exportgraphics(ax3,'figures/fluorescence-1988-2021_lagrangianView.png');

%% Histograms at Depth

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 35 20]);
nbins = 25; 

% Eulerian
[histFreq,histXout] = hist(chloro2D(6,:),nbins);
[histFreq36,histXout36] = hist(chloro2D(36,:),nbins);
[histFreq62,histXout62] = hist(chloro2D(62,:),nbins);

% Lagrangian
[histFreqL,histXoutL] = hist(chloro_lang(2,:),nbins);
[histFreq36L,histXout36L] = hist(chloro_lang(26,:),nbins);
[histFreq62L,histXout62L] = hist(chloro_lang(45,:),nbins);

ax4 = figure;
sgtitle('Histograms of Concentration at Different Depths','FontSize',20);
subplot(2,2,1)
bar(histXout, histFreq/sum(histFreq),'DisplayName','10 db');
hold on
bar(histXout36,histFreq36/sum(histFreq36),'DisplayName','70 db');
bar(histXout62,histFreq62/sum(histFreq62),'DisplayName','122 db');
legend();
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Concentration (Eulerian)');

subplot(2,2,2)
bar(log(histXout), histFreq/sum(histFreq),'DisplayName','10 db');
hold on
bar(log(histXout36),histFreq36/sum(histFreq36),'DisplayName','70 db');
bar(log(histXout62),histFreq62/sum(histFreq62),'DisplayName','122 db');
legend();
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Log-Concentration (Eulerian)');

subplot(2,2,3)
bar(histXoutL, histFreqL/sum(histFreqL),'DisplayName','-100 db');
hold on
bar(histXout36L,histFreq36L/sum(histFreq36L),'DisplayName','-50 db');
bar(histXout62L,histFreq62L/sum(histFreq62L),'DisplayName','-12 db');
legend();
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Log-Concentration (Lagrangian)');

subplot(2,2,4)
bar(log(histXoutL), histFreqL/sum(histFreqL),'DisplayName','-100 db');
hold on
bar(log(histXout36L),histFreq36L/sum(histFreq36L),'DisplayName','-50 db');
bar(log(histXout62L),histFreq62L/sum(histFreq62L),'DisplayName','-12 db');
legend();
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Log-Concentration (Lagrangian)');

exportgraphics(ax4,'figures/hist_chloropig_selectDepths_1989-2021.png');