% ALOHA (A Long-Term Oligotrophic Habitat Assessment; 22° 45'N, 158° 00'W)
% 22.75 N, 158 W

% File to open and examine HOTS T, S, Nitrate Data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',16);

%% Open data file and extract variables
data = importdata('data/hots-T-S-nit.txt').data;

crn = data(:,1);
day = data(:,2);
pressure = data(:,3);
T = data(:,4);
S = data(:,5);
nitrate = data(:,6);

T(T==-9) = NaN;
S(S==-9) = NaN;
nitrate(nitrate==-9) = NaN;

%% Calculate Conservative Temperature, Absolute Salinity, Density

lat = 22.75;
SA = gsw_SA_from_SP(S,pressure,-158,lat);
CT = gsw_CT_from_t(SA,T,pressure);
sigma0 = gsw_sigma0(SA,CT);

%% Reshape CT, SA, sigma0; create Hovmoeller Coordinates

CT2D = reshape(CT,101,[]);
SA2D = reshape(SA,101,[]);
sigma2D = reshape(sigma0,101,[]);

days = reshape(day + datetime(1988,09,30),101,[]);
pres = reshape(pressure,101,[]);

[t_grid,p_grid] = meshgrid(datenum(days(1,:)),pres(:,1));
time = datetime(t_grid(1,:),'ConvertFrom','datenum'); % for non-contour plots

% instead of using p, use sigma0 as vertical coordinate
[~,sigma_grid] = meshgrid(datenum(days(1,:)),sigma2D(:,1));

%% Remove Seasonality in CT, SA, sigma0

CT_rm = movmean(CT2D,10,2,'omitnan');
SA_rm = movmean(SA2D,10,2,'omitnan');
sigma_rm = movmean(sigma2D,10,2,'omitnan');

%% Examine First Days of Data for CT, SA, and sigma

ax1 = figure;
sgtitle('CTD: Temperature & Salinity vs Pressure (HOTS 88/89)');
subplot(1,3,1)
plot(CT(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
hold on
plot(CT(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
plot(CT(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
legend('Location','best');
ylabel('pressure [db]');
xlabel('Conservative Temperature \Theta (Celsius)');

subplot(1,3,2)
plot(SA(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
hold on
plot(SA(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
plot(SA(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
hold off
legend('Location','best');
ylabel('pressure [db]');
xlabel('Absolute Salinity S_A (g/kg)');

subplot(1,3,3)
plot(sigma0(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
hold on
plot(sigma0(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
plot(sigma0(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
hold off
legend('Location','best');
ylabel('pressure [db]');
xlabel('Potential Density Anomaly \sigma^\Theta (kgm^{-3})');

exportgraphics(ax1,'figures/ctd-day-1-SA-CT-sigma.png');

%% SA-CT Plot

p_ref = 1000;
isopycs = 24:0.5:33;

% Adjust size because GSW preset gives a square SA-CT plot
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [16 5 20 20]);

ax2 = figure;
gsw_SA_CT_plot(SA,CT,p_ref,isopycs,'\it{S}\rm_A - {\Theta} diagram');
exportgraphics(ax2,'figures/ctd-SA-CT.png');

%% Conservative Temperature Pressure- and Time-Series (1988-2022): Eulerian View

% New default for 2x1
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [9 2 32 24]);

% Definitions of Mixed Layer Depth (MLD)
% Birol Kara et al. (2000) suggest \Delta T = 0.8 K
% de Boyer Montegut et al (2004) suggest \Delta T = 0.2 C or \Delta \sigma
% = 0.03 kg/m3
% Dong et al (2008) used the above method, so let's try that.
% MLD should be defined with respect to a subsurface layer (Barone, 
% Groeskamp). Thompson et al (2016) use 10m. (10m = level 6)

MLDt = [];
CTTEMP = fillmissing(CT2D(6,:),"previous");
for j = 1:329
    MLDt(j) = find(CT2D(:,j) < CTTEMP(j) - 0.2,1);
end

nb=100

ax3 = figure;
subplot(2,1,1)
h=contourf(t_grid,p_grid,CT2D,linspace(16,28,nb),'LineColor','auto');
hold on
plot(t_grid(1,:),MLDt,'LineWidth',1.5,'Color',[0 0 0]);
set(gca,'Ydir','reverse')
datetickzoom('x','dd/mm/yyyy','keeplimits');
d = colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Conservative Temperature (C)';
xlabel('Time');
ylabel('Pressure [db]');
title('Conservative Temperature \Theta: 1988 - 2021 (Eulerian)');

subplot(2,1,2)
h=contourf(t_grid,p_grid,CT_rm,linspace(16,28,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','dd/mm/yyyy','keeplimits');
d = colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Conservative Temperature (C)';
xlabel('Time');
ylabel('Pressure [db]');
title('Conservative Temperature \Theta: 1988 - 2021 (Eulerian) [12-mth running mean]');

exportgraphics(ax3,'figures/conservativeTemperature-1988-2021_eulerianView.png');

%% Absolute Salinity Pressure- and Time-Series (1988-2022): Eulerian View

ax4 = figure;

subplot(2,1,1)
contourf(t_grid,p_grid,SA2D,linspace(34.6,35.6,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','dd/mm/yyyy','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Absolute Salinity S_A (g kg^{-1})';
xlabel('Time');
ylabel('Pressure [db]');
title('Absolute Salinity S_A: 1988 - 2021 (Eulerian)');

subplot(2,1,2)
contourf(t_grid,p_grid,SA_rm,linspace(34.6,35.6,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','dd/mm/yyyy','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Absolute Salinity S_A (g kg^{-1})';
xlabel('Time');
ylabel('Pressure [db]');
title('Absolute Salinity S_A: 1988 - 2021 (Eulerian) [12-mth running mean]');

exportgraphics(ax4,'figures/absoluteSalinity-1988-2021_eulerianView.png');

%% Density Pressure- and Time-Series (1988-2022): Eulerian View

% Definitions of Mixed Layer Depth (MLD)
% Birol Kara et al. (2000) suggest \Delta T = 0.8 K
% de Boyer Montegut et al (2004) suggest \Delta T = 0.2 C or \Delta \sigma
% = 0.03 kg/m3
% Dong et al (2008) used the above method, so let's try that.

MLD = [];
sigmaTEMP = fillmissing(sigma2D(6,:),"previous");
sigmaTEMP_rm = sigma_rm(6,:);
for j = 1:329
    MLD(j) = find(sigma2D(6:end,j) > sigmaTEMP(j) + 0.03,1);
    MLD_rm(j) = find(sigma_rm(6:end,j) > sigmaTEMP_rm(j) + 0.03,1);
end
% Compensate for starting lower
MLD = MLD + 5;
MLD_rm = MLD_rm + 5;

% Very broad trend
MLD_10yr = movmean(MLD_rm,100);

ax5 = figure;

subplot(2,1,1)
contourf(t_grid,p_grid,sigma2D,linspace(22,25.5,nb),'LineColor','auto');
hold on
plot(t_grid(1,:),MLD,'LineWidth',1.5,'Color',[0 0 0]);
hold off
set(gca,'Ydir','reverse')
datetickzoom('x','dd/mm/yyyy','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'Potential Density Anomaly \sigma^\Theta (kg m^{-3})';
xlabel('Time');
ylabel('Pressure [db]');
title('Potential Density Anomaly \sigma^\Theta: 1988 - 2021 (Eulerian)');

subplot(2,1,2)
contourf(t_grid,p_grid,sigma_rm,linspace(22,25.5,nb),'LineColor','auto');
hold on
plot(t_grid(1,:),MLD_rm,'LineWidth',1.5,'Color',[0 0 0]);
plot(t_grid(1,:),MLD_10yr,'LineWidth',1.5,'Color',[0.3 0.3 0.3]);
hold off
set(gca,'Ydir','reverse')
datetickzoom('x','dd/mm/yyyy','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'Potential Density Anomaly \sigma^\Theta (kg m^{-3})';
xlabel('Time');
ylabel('Pressure [db]');
title('Potential Density Anomaly \sigma^\Theta: 1988 - 2021 (Eulerian) [12-mth running mean]');

exportgraphics(ax5,'figures/density-1988-2021_eulerianView.png');

%% Buoyancy Time- and Depth-Series: Eulerian

% Return to default view for Pressure- and time-series
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);

[N2,pmid] = gsw_Nsquared(SA2D,CT2D,pres,lat);
[N2_rm,pmid_rm] = gsw_Nsquared(SA_rm,CT_rm,pres,lat);

[tmid_grid,pmid_grid] = meshgrid(datenum(days(1,:)),pmid(:,1));

ax6 = figure;

subplot(2,1,1)
contourf(tmid_grid,pmid_grid,real(log10(N2)),linspace(-8.5,-3,nb),'LineColor','auto');
hold on
plot(tmid_grid(1,:),MLD,'LineWidth',1.5,'Color',[0 0 0]);
hold off
set(gca,'Ydir','reverse')
datetickzoom('x','dd/mm/yyyy','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Buoyancy N^2 (s^{-2})';
xlabel('Time');
ylabel('Pressure [db]');
title('Buoyancy N^2: 1988 - 2021 (Eulerian)');

subplot(2,1,2)
contourf(tmid_grid,pmid_grid,real(log10(N2_rm)),linspace(-8.5,-3,nb),'LineColor','auto');
hold on
plot(tmid_grid(1,:),MLD_rm,'LineWidth',1.5,'Color',[0 0 0]);
hold off
set(gca,'Ydir','reverse')
datetickzoom('x','dd/mm/yyyy','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Buoyancy N^2 (s^{-2})';
xlabel('Time');
ylabel('Pressure [db]');
title('Buoyancy N^2: 1988 - 2021 (Eulerian) [12-mth running mean]');

exportgraphics(ax6,'figures/buoyancy-1988-2021_eulerianView.png');

%% FFT on N2

% Fill N2, calculate FFT
for i = 1:100
    N2filled(i,:) = fillmissing(N2(i,:),'previous');
    y(i,:) = nufft(N2filled(i,:),tmid_grid(i,:));
end

meanN2fft = mean(y);

n = length(time);
f = (0:n-1)/n;

ax7 = figure;
plot(f,abs(y),'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hold on
plot(f,abs(meanN2fft),'Color','red','DisplayName','Mean N^2 Across Depth');
hold off
legend();
xlabel('Normalised Frequency');
ylabel('Power Density (N^4) [s^{-4}]');
title('FFT: N^2(p,t) at HOTS (1989-2021)');
exportgraphics(ax7,'figures/buoyancy-1988-2021_fft.png');

%% FFT on MLD

% Fill N2, calculate FFT
y_mld = nufft(MLD,tmid_grid(i,:));

ax8 = figure;
yyaxis left
plot(f,abs(y_mld),'Color',[0.8 0.8 0.8],'DisplayName','Mixed Layer Depth');
xlabel('Normalised Frequency');
ylabel('MLD: Power Density (p^2) [Pa^2]');
hold on
yyaxis right
plot(f,abs(meanN2fft),'Color','red','DisplayName','Mean N^2 Across Depth');
hold off
legend();
xlabel('Normalised Frequency');
ylabel('N^2: Power Density (N^4) [s^{-4}]');
title('FFT of MLD and depth-averaged N^2 (HOTS, 1989-2021)');
exportgraphics(ax8,'figures/mld-1988-2021_fft.png');

%% MLD and depth-averaged N^2 are correlated

% Check this.

%% Density vs chl

load("datafiles\chloro.mat");
ax9 = figure;

subplot(2,1,1)
contourf(t_grid,-sigma_grid,chloro2D,linspace(0,1.4,nb),'LineColor','auto');
datetickzoom('x','dd/mm/yyyy','keeplimits');
d = colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Chloropigment (ug/L)';
xlabel('Time');
ylabel('Potential Density Anomaly $\sigma_0$');
title('Chloropigment (ug/L): 1988 - 2021 (Eulerian, Sigma Coordinates)');

subplot(2,1,2)
contourf(t_grid,-p_grid,chloro2D,linspace(0,1.4,nb),'LineColor','auto');
datetickzoom('x','dd/mm/yyyy','keeplimits');
d = colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'Chloropigment (ug/L)';
xlabel('Time');
ylabel('Pressure [db]');
title('Chloropigment (ug/L): 1988 - 2021 (Eulerian, Pressure Coordinates)');

exportgraphics(ax8,'figures/chl_sigmaCoords.png');

%% SA-CT with Chl-a: discretize Chl-a data to map onto 100 colours from Hovmoeller Diagram
testDisc = discretize(chloro2D_n,100);
testDisc(isnan(testDisc)) = 101;
d = [d;0 0 0];

%% SA-CT with Chl-a: set up density coordinates

SA_sig = linspace(min(min(SA2D))-0.1,max(max(SA2D))+0.1,nb);
CT_sig = linspace(min(min(CT2D))-1,max(max(CT2D)),nb);
[SA_sigg,CT_sigg] = meshgrid(SA_sig,CT_sig);
sigma_lines = gsw_sigma0(SA_sigg,CT_sigg);

%% SA-CT with Chl-a: Make the Diagram

ax9 = figure;
contour(SA_sigg,CT_sigg,sigma_lines,'k-','ShowText','on','LabelSpacing',800);
hold on
gscatter(SA2D,CT2D,d(testDisc(:,1),:),d,'.',8,'off');
hold off
xlabel('Absolute Salinity S_A [g kg^{-1}]');
ylabel('Conservative Temperature \Theta [C]');
exportgraphics(ax9,'figures/CT_SA_withChl.png');