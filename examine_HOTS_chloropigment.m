% File to open and examine HOTS chloropigment data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',16);

%% Open data file and extract variables; add Barone's routines to path
data = importdata('data/hots-chloropigment.txt').data;

crn = data(:,1);
day = data(:,2);
pressure = data(:,3);
chloro = data(:,4);

% Set -9 = NaN; remove negative values
chloro(chloro==-9) = NaN;
chloro(chloro<0) = 0;

addpath("baroneRoutines\");

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
time = datetime(t_grid(1,:),'ConvertFrom','datenum');

save('datafiles\chloro',"chloro2D","pres","t_grid"',"p_grid","time");


ax2 = figure;
contourf(t_grid,p_grid,chloro2D,0:0.14:1.4,'LineColor','auto');
set(gca,'Ydir','reverse')
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Eulerian)');

exportgraphics(ax2,'figures/fluorescence-1988-2021_eulerianView.png');

%% Normalised Chloropigment Depth- and Time-Series (1988-2020): Eulerian

for j=1:329
    chloro2D_n(:,j) = chloro2D(:,j)/max(chloro2D(:,j));
end
save("datafiles\chloro.mat","chloro2D_n",'-append');

ax2a = figure;
contourf(t_grid,p_grid,chloro2D_n,'LineColor','auto');
set(gca,'Ydir','reverse')
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'chloropigment, normalised relative to DCM at each time';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Eulerian, Normalised)');

exportgraphics(ax2a,'figures/fluorescence_norm-1988-2021_eulerianView.png');

% Test normalisation
assert(max(max(chloro2D_n)) == 1); % Throws error if not equal to one
assert(min(min(chloro2D_n)) == 0); % Throws error if not equal to zero

%% Kurtosis and Skewness across depth for (normalised) chloropigment depth- and time-series (Eulerian)

kurt_chl = kurtosis(chloro2D);
skew_chl = skewness(chloro2D);

kurt_chl_rm = movmean(kurt_chl,10);
skew_chl_rm = movmean(skew_chl,10);
% The normalised data of course has the same kurtosis and skewness
% kurt_chl_n = kurtosis(chloro2D_n);
% skew_chl_n = skewness(chloro2D_n);

ax2b = figure;
plot(time,kurt_chl,'DisplayName','Kurtosis');
hold on
plot(time,kurt_chl_rm,'DisplayName','~12mth running-mean (10-point centred moving average');
% plot(t_grid(1,:),kurt_chl_n,'DisplayName','Kurtosis (Norm)');
yline(3,':','DisplayName','Normal Distribution');
hold off
legend();
title('Kurtosis: Chloropigments, 1988-2021 (Eulerian)');

exportgraphics(ax2b,'figures/fluorescence_norm-1988-2021_eulerianKurtosis.png');

ax2c = figure;
plot(time,skew_chl,'DisplayName','Skewness');
hold on
plot(time,skew_chl_rm,'DisplayName','~12mth running mean (10-point centred moving average)');
hold off
% plot(t_grid(1,:),skew_chl_n,'DisplayName','Skewness (norm)');
legend();
datetick('x','yyyy mmm','keeplimits');
title('Skewness: Chloropigments, 1988-2021 (Eulerian)');

exportgraphics(ax2c,'figures/fluorescence_norm-1988-2021_eulerianSkewness.png');


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
title('Chloropigment: 1988 - 2021 (Lagrangian)');

exportgraphics(ax3,'figures/fluorescence-1988-2021_lagrangianView.png');

%% Lagrangian, Normalised

for j=1:329
    chloro_lang_n(:,j) = chloro_lang(:,j)/max(chloro_lang(:,j));
end

ax3a = figure;
contourf(t_lang_grid,p_lang_grid,chloro_lang_n,'LineColor','auto');
set(gca,'Ydir','reverse')
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'chloropigment, normalised relative to DCM';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Lagrangian, Normalised)');

exportgraphics(ax3a,'figures/fluorescence_norm-1988-2021_lagrangianView.png');

%% Kurtosis and Skewness across depth for normalised chloropigment depth- and time-series (Lagrangian)

kurt_chl_lang = kurtosis(chloro_lang);
kurt_chl_lang_rm = movmean(kurt_chl_lang,10);
skew_chl_lang = skewness(chloro_lang);
skew_chl_lang_rm = movmean(skew_chl_lang,10);

% Again the normalised data exhibits the same kurtosis and skew as the
% non-normalised data
% kurt_chl_lang_n = kurtosis(chloro_lang_n);
% skew_chl_lang_n = skewness(chloro_lang_n);

ax3b = figure;
plot(time,kurt_chl_lang,'DisplayName','Kurtosis');
hold on
plot(time,kurt_chl_lang_rm,'DisplayName','~12mth running mean (10-point centred moving average');
% plot(t_grid(1,:),kurt_chl_lang_n,'DisplayName','Kurtosis (norm)');
yline(3,':','DisplayName','Normal Distribution');
hold off
legend();
datetick('x','yyyy mmm','keeplimits');
title('Kurtosis: Chloropigments, 1988-2021 (Lagrangian)');

exportgraphics(ax3b,'figures/fluorescence_norm-1988-2021_lagrangianKurtosis.png');

ax3c = figure;
plot(time,skew_chl_lang,'DisplayName','Skewness');
hold on
plot(time,skew_chl_lang_rm,'DisplayName','~12mth running mean (10-point centred moving average)');
% plot(t_grid(1,:),skew_chl_lang_n,'DisplayName','Skewness (norm)');
hold off
legend();
datetick('x','yyyy mmm','keeplimits');
title('Skewness: Chloropigments, 1988-2021 (Lagrangian)');

exportgraphics(ax3c,'figures/fluorescence_norm-1988-2021_lagrangianSkewness.png');

%% Histograms at Depth

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 35 20]);
nbins = 25; 

% Eulerian
[histFreq,histXout] = hist(chloro2D(6,:),nbins);
[histFreq36,histXout36] = hist(chloro2D(36,:),nbins);
[histFreq62,histXout62] = hist(chloro2D(62,:),nbins);

% Lagrangian
[histFreqL,histXoutL] = hist(chloro_lang(2,:),nbins);
[histFreq26L,histXout26L] = hist(chloro_lang(26,:),nbins);
[histFreq45L,histXout45L] = hist(chloro_lang(45,:),nbins);
[histFreq12L,histXout12L] = hist(chloro_lang(58,:),nbins);
[histFreq50L,histXout50L] = hist(chloro_lang(77,:),nbins);

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
bar(histXout26L,histFreq26L/sum(histFreq26L),'DisplayName','-50 db');
bar(histXout45L,histFreq45L/sum(histFreq45L),'DisplayName','-12 db');
bar(histXout12L,histFreq12L/sum(histFreq12L),'DisplayName','+12 db');
bar(histXout50L,histFreq50L/sum(histFreq50L),'DisplayName','+50 db');
legend();
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Log-Concentration (Lagrangian)');

subplot(2,2,4)
bar(log(histXoutL), histFreqL/sum(histFreqL),'DisplayName','-100 db');
hold on
bar(log(histXout26L),histFreq26L/sum(histFreq26L),'DisplayName','-50 db');
bar(log(histXout45L),histFreq45L/sum(histFreq45L),'DisplayName','-12 db');
bar(log(histXout12L),histFreq12L/sum(histFreq12L),'DisplayName','+12 db');
bar(log(histXout50L),histFreq50L/sum(histFreq50L),'DisplayName','+50 db');
legend();
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Log-Concentration (Lagrangian)');

exportgraphics(ax4,'figures/hist_chloropig_selectDepths_1989-2021.png');