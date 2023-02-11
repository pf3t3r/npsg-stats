% File to open and examine HOTS chloropigment data (1988 - 2021).

% Clear/close unnecessary code, variables, etc.
close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 28 12]);
set(0,'defaultAxesFontSize',12);

%% Open data file and extract variables; add Barone's routines to path

% Fluorescence Data
data = importdata('data/hots-chloropigment.txt').data;
data1 = importdata('data\hots-chl-T-S-nit-1000.txt').data;
data256 = importdata('data\hots-chl-T-S-nit-256.txt').data;
d_300 = importdata('data\hots-T-S-chl-nit-300.txt').data;

% Bottle Data
cpig_bot = importdata('data\hot-bottle.txt').data;

%% Check CTD Data for Completeness 
% Do they go full depth? Check that for each cast the pressure goes down to the specified value

% 256 db
noOfZeroDbs256 = find(~data256(:,3));
assert(length(noOfZeroDbs256)==329);
noOf256Dbs256 = find(data256(:,3) == 256);
assert(length(noOf256Dbs256)==329);

% 300 db: this test will FAIL so it is commented out
% noOfZeroDbs300 = find(~d_300(:,3));
% assert(length(noOfZeroDbs300)==329);
% noOf300Dbs300 = find(d_300(:,3) == 300);
% assert(length(noOf300Dbs300)==329);

% There are two CTDs which do not go down to 300 db
% Therefore we must remove them
idIncompleteCTDpres = find(diff(d_300(:,3)) < 2 & diff(d_300(:,3)) > -300);
idStartRemoval = idIncompleteCTDpres - d_300(idIncompleteCTDpres,3)/2;

% The two incomplete casts follow each other
d_300 = d_300([1:idStartRemoval(1)-1 idIncompleteCTDpres(2)+1:end],:);

% Test again: Two casts removed => length = 327
% This test should be passed
noOfZeroDbs300 = find(~d_300(:,3));
assert(length(noOfZeroDbs300)==327);
noOf300Dbs300 = find(d_300(:,3) == 300);
assert(length(noOf300Dbs300)==327);

%% Extract CTD Data, Assign to Variables

% 200db
crn = data(:,1);
day = data(:,2);
pressure = data(:,3);
chloro = data(:,4);
chloro(chloro==-9) = NaN; % Set -9 = NaN;
chloro(chloro<0) = 0; % set negative values = 0

% 256db
crn2 = data256(:,1);
day2 = data256(:,2);
pressure2 = data256(:,3);
chloro2 = data256(:,6);
chloro2(chloro2==-9) = NaN; % Set -9 = NaN;
chloro2(chloro2<0) = 0; % set negative values = 0

% 300 db
crn3 = d_300(:,1);
day3 = d_300(:,2);
pres3 = d_300(:,3);
chloro3 = d_300(:,6);
chloro3(chloro3==-9) = NaN; % Set -9 = NaN;
chloro3(chloro3<0) = 0; % set negative values = 0

% 1000db
crn1 = data1(:,1);
day1 = data1(:,2);
pressure1 = data1(:,3);
chloro1 = data1(:,6);
chloro1(chloro1==-9) = NaN; % Set -9 = NaN;
chloro1(chloro1<0) = 0; % remove negative values

addpath("baroneRoutines\");

nb = 100;

%% Extract Bottle Data

botPressure = cpig_bot(:,4);
botChl = cpig_bot(:,5);

% Grid Data


%% Examine First Days of Data
% Not really required. Can delete later.

% ax1 = figure;
% subplot(1,3,1)
% plot(chloro(1:101),-pressure(1:101),'DisplayName','CRN 1: Nov 88');
% hold on
% plot(chloro(102:202),-pressure(102:202),'DisplayName','CRN 2: Dec 88');
% plot(chloro(203:303),-pressure(203:303),'DisplayName','CRN 3: Jan 89');
% hold off
% legend('Location','best');
% ylabel('pressure [db]');
% xlabel('chl [ug/L]');
% title('chl(p) [200db]');
% 
% subplot(1,3,2)
% plot(chloro1(1:501),-pressure1(1:501),'DisplayName','CRN 1: Nov 88');
% hold on
% plot(chloro1(502:1002),-pressure1(502:1002),'DisplayName','CRN 2: Dec 88');
% plot(chloro1(1003:1503),-pressure1(1003:1503),'DisplayName','CRN 3: Jan 89');
% hold off
% legend('Location','best');
% ylabel('pressure [db]');
% xlabel('chl [ug/L]');
% title('chl(p) [1000db]');
% 
% subplot(1,3,3)
% plot(chloro2(1:129),-pressure2(1:129),'DisplayName','CRN 1: Nov 88');
% hold on
% plot(chloro2(130:258),-pressure2(130:258),'DisplayName','CRN 2: Dec 88');
% plot(chloro2(259:387),-pressure2(259:387),'DisplayName','CRN 3: Jan 89');
% hold off
% legend('Location','best');
% ylabel('pressure [db]');
% xlabel('chl [ug/L]');
% title('chl(p) [256db]');
% 
% exportgraphics(ax1,'figures/ctd-day-1.png');

%% Reshape into 2D Matrix: 200db

chloro200 = reshape(chloro,101,[]);
days200 = day + datetime(1988,09,30);
days200 = reshape(days200,101,[]);
pres200 = reshape(pressure,101,[]);

save('datafiles\chloro',"chloro200",'days200',"pres200");

%% Reshape into 2D Matrix: 256db

chloro256 = reshape(chloro2,129,[]);
days256 = day2 + datetime(1988,09,30);
days256 = reshape(days256,129,[]);
pres256 = reshape(pressure2,129,[]);

save('datafiles\chloro',"chloro256",'days256',"pres256",'-append');

%% Reshape into 2D Matrix: 300 db

chloro300 = reshape(chloro3,151,[]);
days300 = day3 + datetime(1988,09,30);
days300 = reshape(days300,151,[]);
pres300 = reshape(pres3,151,[]);

save('datafiles\chloro',"chloro300",'days300',"pres300",'-append');

%% Reshape into 2D Matrix: 1000db

% There are two casts where the CTD does not go to the full depth. I would
% like to fill these missing values with NaNs and be able to plot a full
% depth Hovmoeller diagram of chl. But I need to return to this later.

% % CRN = 329
% noOfZeroDbs = find(~pressure1);
% assert(length(noOfZeroDbs)==329);
% 
% counter = 0;
% counterPression = [];
% for i=1:length(pressure1)-1
%     if diff(pressure1(i:i+1)) < 2 && diff(pressure1(i:i+1)) > -999
%         counter = counter + 1;
%         counterPression = [counterPression i];
%         disp(i);
%     end
% end
% 
% % Insert NaNs
% copyPressure1 = pressure1;
% Nangroup1 = NaN((1000-266)/2,1);
% Nangroup2 = NaN((1000-270)/2,1);
% copyPressure1(counterPression(1))
% 
% % copyPressure1 = [copyPressure1(1:counterPression(1)),Nangroup1,copyPressure1(counterPression(1)+1:end)];
% 
% copytest = [copyPressure1(1:counterPression(1))',Nangroup1'];
% copytest = [copytest,copyPressure1(counterPression(1)+1:counterPression(2))'];
% copytest = [copytest,Nangroup2'];
% copytest = [copytest,copyPressure1(counterPression(2)+1:end)'];
% %     copyPressure1(counterPression(1)+1:counterPression(2)),Nangroup2,...
% %     copyPressure1(counterPression(2)+1:end)];
% 
% % chloro1_2D = reshape(chloro1,501,[]);
% % days1 = day1 + datetime(1988,09,30);
% % days1 = reshape(days1,501,[]);
% pres1 = reshape(copytest,501,[]);

%% Apply Gridding.

[tgrid200,pgrid200] = meshgrid(datenum(days200(1,:)),pres200(:,1));
time200 = datetime(tgrid200(1,:),'ConvertFrom','datenum');

[tgrid256,pgrid256] = meshgrid(datenum(days256(1,:)),pres256(:,1));
time256 = datetime(tgrid256(1,:),'ConvertFrom','datenum');

[tgrid300,pgrid300] = meshgrid(datenum(days300(1,:)),pres300(:,1));
time300 = datetime(tgrid300(1,:),'ConvertFrom','datenum');

save('datafiles\chloro',...
    "tgrid200"',"pgrid200","time200",...
    "tgrid256","pgrid256","time256",...
    "tgrid300","pgrid300","time300",'-append');

%% Chloropigment: Eulerian Hovmoeller (1988-2022) [200, 256,  db]

% 200 db
ax2 = figure;
contourf(tgrid200,pgrid200,chloro200,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: Eulerian (1988-2021) [200 db]');

exportgraphics(ax2,'figures/fluorescence-1988-2021_eulerian.png');

% 256 db
ax3 = figure;
contourf(tgrid256,pgrid256,chloro256,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: Eulerian (1988-2021) [256db]');

exportgraphics(ax3,'figures/fluorescence-1988-2021_256_eulerian.png');

% 300 db
ax3a = figure;
contourf(tgrid300,pgrid300,chloro300,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: Eulerian (1988 - 2021) [300 db]');

exportgraphics(ax3a,'figures/fluorescence-1988-2021_300_eulerian.png');

clear ax2 ax3 ax3a;
%% Chloropigment (Normalised), Eulerian Hovmoeller (1988-2020)

chloro200n = zeros(101,329);
chloro256n = zeros(129,329);
chloro300n = zeros(151,327);

for j=1:329
    chloro200n(:,j) = chloro200(:,j)/max(chloro200(:,j));
    chloro256n(:,j) = chloro256(:,j)/max(chloro256(:,j));
end
for k=1:327
    chloro300n(:,k) = chloro300(:,k)/max(chloro300(:,k));
end

save("datafiles\chloro.mat","chloro200n",'chloro256n','chloro300n','-append');
clear j k;

% Test that values were normalised properly
assert(max(max(chloro200n)) == 1); % Throws error if not equal to one
assert(min(min(chloro200n)) == 0); % Throws error if not equal to zero
assert(max(max(chloro256n)) == 1); % Throws error if not equal to one
assert(min(min(chloro256n)) == 0); % Throws error if not equal to zero
assert(max(max(chloro300n)) == 1); % Throws error if not equal to one
assert(min(min(chloro300n)) == 0); % Throws error if not equal to zero

ax4 = figure;
contourf(tgrid200,pgrid200,chloro200n,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment, normalised relative to DCM at each time';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Eulerian, Normalised)');

exportgraphics(ax4,'figures/fluorescence_norm-1988-2021_eulerianView.png');

ax4a = figure;
contourf(tgrid256,pgrid256,chloro256n,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment, normalised relative to DCM at each time';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Eulerian, Normalised) [256 db]');

exportgraphics(ax4a,'figures/fluorescence_norm-1988-2021_eulerian256.png');

ax4b = figure;
contourf(tgrid300,pgrid300,chloro300n,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment, normalised relative to DCM at each time';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Eulerian, Normalised) [300 db]');
exportgraphics(ax4b,'figures/fluorescence_norm-1988-2021_eulerian300.png');

clear ax4 ax4a ax4b;
%% Kurtosis and Skewness across depth for (normalised) chloropigment depth- and time-series (Eulerian)
% 
% kurt_chl = kurtosis(chloro2D);
% skew_chl = skewness(chloro2D);
% 
% kurt_chl_rm = movmean(kurt_chl,10,'omitnan');
% skew_chl_rm = movmean(skew_chl,10,'omitnan');
% % The normalised data of course has the same kurtosis and skewness
% % kurt_chl_n = kurtosis(chloro2D_n);
% % skew_chl_n = skewness(chloro2D_n);
% 
% ax5 = figure;
% plot(time,kurt_chl,'DisplayName','Kurtosis');
% hold on
% plot(time,kurt_chl_rm,'DisplayName','~12mth running-mean (10-point centred moving average');
% % plot(t_grid(1,:),kurt_chl_n,'DisplayName','Kurtosis (Norm)');
% yline(3,':','DisplayName','Normal Distribution');
% hold off
% legend();
% title('Kurtosis: Chloropigments, 1988-2021 (Eulerian)');
% 
% exportgraphics(ax5,'figures/fluorescence_norm-1988-2021_eulerianKurtosis.png');
% 
% ax6 = figure;
% plot(time,skew_chl,'DisplayName','Skewness');
% hold on
% plot(time,skew_chl_rm,'DisplayName','~12mth running mean (10-point centred moving average)');
% hold off
% % plot(t_grid(1,:),skew_chl_n,'DisplayName','Skewness (norm)');
% legend();
% datetickzoom('x','yyyy mmm','keeplimits');
% title('Skewness: Chloropigments, 1988-2021 (Eulerian)');
% 
% exportgraphics(ax6,'figures/fluorescence_norm-1988-2021_eulerianSkewness.png');


%% Reformulate Eulerian Data into Lagrangian View

val200 = zeros(1,329); idx200 = zeros(1,329); 
val256 = zeros(1,329); idx256 = zeros(1,329);
val300 = zeros(1,327); idx300 = zeros(1,327);

% Find the DCM
for i=1:329
    [val200(i),idx200(i)] = max(chloro200(:,i)); % 200 db
    [val256(i),idx256(i)] = max(chloro256(:,i));
end
for j=1:327
    [val300(j),idx300(j)] = max(chloro300(:,j)); % 300 db
end
clear i j;

pres200l = zeros(101,329);
pres256l = zeros(129,329);
pres300l = zeros(151,327);

% Put pressure in terms of DCM
for i = 1:329
    pres200l(:,i) = pres200(:,i) - pres200(idx200(i),i);
    pres256l(:,i) = pres256(:,i) - pres256(idx256(i),i);
end
for j = 1:327
    pres300l(:,j) = pres300(:,j) - pres300(idx300(j),j);
end
clear i j;

save("datafiles\chloro.mat",'pres200l',"pres256l","pres300l",'-append');

%%
% Put fluorescence in terms of DCM
% Shift the original fluorescence data such that the DCM is centred
chloro200l = zeros(101,329);
midpt200 = 51;
offset200 = midpt200 - idx200;

chloro256l = zeros(129,329);
midpt256 = 65;
offset256 = midpt256 - idx256;

chloro300l = zeros(151,327);
midpt300 = 76;
offset300 = midpt300 - idx300;

% 200 db
for i = 1:329
    chloro200l(:,i) = circshift(chloro200(:,i),offset200(i));
    if offset200(i) > -1 && offset200(i) < 40
        disp(i);
        chloro200l(1:offset200(i),i) = NaN;
    elseif offset200(i) == -1
        chloro200l(end,i) = NaN;
    elseif offset200(i) < -1 && offset200(i) > -40
        disp(i);
        chloro200l((end+offset200(i)):end,i) = NaN;
    elseif abs(offset200(i)) > 40
        chloro200l(:,i) = NaN;
    end
end

% 256 db
for i = 1:329
    chloro256l(:,i) = circshift(chloro256(:,i),offset256(i));
    if offset256(i) > -1 && offset256(i) < 40
        disp(i);
        chloro256l(1:offset256(i),i) = NaN;
    elseif offset256(i) == -1
        chloro256l(end,i) = NaN;
    elseif offset256(i) < -1 && offset256(i) > -40
        disp(i);
        chloro256l((end+offset256(i)):end,i) = NaN;
    elseif abs(offset256(i)) > 40
        chloro256l(:,i) = NaN;
    end
end

% 300 db
for i = 1:327
    chloro300l(:,i) = circshift(chloro300(:,i),offset300(i));
    if offset300(i) > -1 && offset300(i) < 40
        disp(i);
        chloro300l(1:offset300(i),i) = NaN;
    elseif offset300(i) == -1
        chloro300l(end,i) = NaN;
    elseif offset300(i) < -1 && offset300(i) > -40
        disp(i);
        chloro300l((end+offset300(i)):end,i) = NaN;
    elseif abs(offset300(i)) > 40
        chloro300l(:,i) = NaN;
    end
end

save("datafiles\chloro.mat",'chloro200l',"chloro256l","chloro300l",'-append');


%%
% save('datafiles\chloro', ...
%     "chloro200","pres200","t_grid"',"pgrid_200","time200", ...
%     'chloro256n','chloro256',...
%     "tgrid_256","pgrid_256","time256");

%%
% Create meshgrid for time and pressure in Lagrangian view
[tgrid200l,pgrid200l] = meshgrid(datenum(days200(1,:)),pres200(:,1)-100);
[tgrid256l,pgrid256l] = meshgrid(datenum(days256(1,:)),pres256(:,1)-128);
[tgrid300l,pgrid300l] = meshgrid(datenum(days300(1,:)),pres300(:,1)-150);

save('datafiles\chloro',...
    'tgrid200l','pgrid200l',...
    'tgrid256l','pgrid256l',...
    'tgrid300l','pgrid300l',...
    'pgrid256l','-append');

%% Chloropigment: Lagrangian Hovmoeller (1988-2022)

% 200 db
ax7 = figure;
contourf(tgrid200l,pgrid200l,chloro200l,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Lagrangian)');
exportgraphics(ax7,'figures/fluorescence-1988-2021_lagrangianView.png');

% 256 db
ax7a = figure;
contourf(tgrid256l,pgrid256l,chloro256l,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
ylim([-110 110]);
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Lagrangian) [256 db]');
exportgraphics(ax7a,'figures/fluorescence-1988-2021_lagrangian256.png');

% 300 db
ax7b = figure;
contourf(tgrid300l,pgrid300l,chloro300l,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment (fluorescence) [ug/L]';
ylim([-110 110]);
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Lagrangian) [300 db]');
exportgraphics(ax7b,'figures/fluorescence-1988-2021_lagrangian300.png');

clear ax7 ax7a ax7b;
%% Chloropigment (Normalised): Lagrangian Hovmoeller (1988-2022)

chloro200nl = zeros(101,329);
chloro256nl = zeros(129,329);
chloro300nl = zeros(151,327);

for j=1:329
    chloro200nl(:,j) = chloro200l(:,j)/max(chloro200l(:,j));
    chloro256nl(:,j) = chloro256l(:,j)/max(chloro256l(:,j));
end
for j=1:327
    chloro300nl(:,j) = chloro300l(:,j)/max(chloro300l(:,j));
end
clear j;

save('datafiles\chloro',...
    'chloro200nl','chloro256nl','chloro300nl', ...
    '-append');

% 200 db
ax8 = figure;
contourf(tgrid200l,pgrid200l,chloro200nl,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment, normalised relative to DCM';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Lagrangian, Normalised)');
exportgraphics(ax8,'figures/fluorescence_norm-1988-2021_lagrangianView.png');

% 256 db
ax8a = figure;
contourf(tgrid256l,pgrid256l,chloro256nl,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment, normalised relative to DCM';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Lagrangian, Normalised, 256 db)');
exportgraphics(ax8a,'figures/fluorescence_norm-1988-2021_lagrangian256.png');

% 300 db
ax8b = figure;
contourf(tgrid300l,pgrid300l,chloro300nl,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment, normalised relative to DCM';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment (Normalised): Lagrangian (1988 - 2021) [300 db]');
exportgraphics(ax8b,'figures/fluorescence_norm-1988-2021_lagrangian300.png');

clear ax8 ax8a ax8b;
%% Kurtosis and Skewness across depth for normalised chloropigment depth- and time-series (Lagrangian)
% Remove later.

% 
% kurt_chl_lang = kurtosis(chloro_lang);
% kurt_chl_lang_rm = movmean(kurt_chl_lang,10,'omitnan');
% skew_chl_lang = skewness(chloro_lang);
% skew_chl_lang_rm = movmean(skew_chl_lang,10,'omitnan');
% 
% % Again the normalised data exhibits the same kurtosis and skew as the
% % non-normalised data
% % kurt_chl_lang_n = kurtosis(chloro_lang_n);
% % skew_chl_lang_n = skewness(chloro_lang_n);
% 
% ax9 = figure;
% plot(time,kurt_chl_lang,'DisplayName','Kurtosis');
% hold on
% plot(time,kurt_chl_lang_rm,'DisplayName','~12mth running mean (10-point centred moving average');
% % plot(t_grid(1,:),kurt_chl_lang_n,'DisplayName','Kurtosis (norm)');
% yline(3,':','DisplayName','Normal Distribution');
% hold off
% legend();
% datetickzoom('x','yyyy mmm','keeplimits');
% title('Kurtosis: Chloropigments, 1988-2021 (Lagrangian)');
% 
% exportgraphics(ax9,'figures/fluorescence_norm-1988-2021_lagrangianKurtosis.png');
% 
% ax10 = figure;
% plot(time,skew_chl_lang,'DisplayName','Skewness');
% hold on
% plot(time,skew_chl_lang_rm,'DisplayName','~12mth running mean (10-point centred moving average)');
% % plot(t_grid(1,:),skew_chl_lang_n,'DisplayName','Skewness (norm)');
% hold off
% legend();
% datetickzoom('x','yyyy mmm','keeplimits');
% title('Skewness: Chloropigments, 1988-2021 (Lagrangian)');
% 
% exportgraphics(ax10,'figures/fluorescence_norm-1988-2021_lagrangianSkewness.png');

%% Histograms at Depth

nbins = 25; 

ax11 = figure;
subplot(2,2,1)
histogram(chloro200(6,:),nbins,'DisplayName','10 db');
hold on
histogram(chloro200(36,:),nbins,'DisplayName','70 db');
histogram(chloro200(62,:),nbins,'DisplayName','122 db');
hold off
legend('Location','best');
xlim([0 1.65]);
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Concentration (Eulerian)');

subplot(2,2,2)
histogram(log(chloro200(6,:)),nbins,'DisplayName','10 db');
hold on
histogram(log(chloro200(36,:)),nbins,'DisplayName','70 db');
histogram(log(chloro200(62,:)),nbins,'DisplayName','122 db');
hold off
legend('Location','best');
xlim([-9 0.2]);
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Log-Concentration (Eulerian)');

subplot(2,2,3)
histogram(chloro200l(2,:),nbins,'DisplayName','-100 db');
hold on
histogram(chloro200l(26,:),nbins,'DisplayName','-50 db');
histogram(chloro200l(45,:),nbins,'DisplayName','-12 db');
histogram(chloro200l(58,:),nbins,'DisplayName','+12 db');
histogram(chloro200l(77,:),nbins,'DisplayName','+50 db');
hold off
legend('Location','best');
xlim([0 1.65]);
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Concentration (Lagrangian)');

subplot(2,2,4)
histogram(log(chloro200l(2,:)),nbins,'DisplayName','-100 db');
hold on
histogram(log(chloro200l(26,:)),nbins,'DisplayName','-50 db');
histogram(log(chloro200l(45,:)),nbins,'DisplayName','-12 db');
histogram(log(chloro200l(58,:)),nbins,'DisplayName','+12 db');
histogram(log(chloro200l(77,:)),nbins,'DisplayName','+50 db');
hold off
legend('Location','best');
xlim([-9 0.2]);
xlabel('chloropigment (\mu g L^{-1})');
ylabel('Frequency');
title('Frequency of Log-Concentration (Lagrangian)');

exportgraphics(ax11,'figures/hist_chloropig_selectDepths_1989-2021.png');
clear ax11;

%% DCM Hovmoeller: Seasonality Removed

chloro256_rm = movmean(chloro256,10,2,'omitnan');
chloro256_lang_rm = movmean(chloro256l,10,2,'omitnan');

ax12 = figure;
contourf(tgrid256l,pgrid256l,chloro256_lang_rm,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment, seasonality removed';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Lagrangian, Seasonality Removed, 256 db)');

exportgraphics(ax12,'figures/fluorescence_1988-2021_lagrangian256_seasonalityRemoved.png');

ax13 = figure;
contourf(tgrid256l,pgrid256l,chloro256_rm,linspace(0,1.4,nb),'LineColor','auto');
set(gca,'Ydir','reverse')
datetickzoom('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',nb)));
c = colorbar;
c.Label.String = 'chloropigment, seasonality removed';
xlabel('Time');
ylabel('Depth [db]');
title('Chloropigment: 1988 - 2021 (Eulerian, Seasonality Removed, 256 db)');

exportgraphics(ax13,'figures/fluorescence_1988-2021_eulerian256_seasonalityRemoved.png');
clear ax12 ax13;

%% FFT on DCM

% For now, run the weird FFT on the seasonally-averaged data. Maybe this
% gives something interesting.

y = zeros(129,329);

% Maybe filling missing data not needed here.
for i = 1:129
    %chlorofilled(i,:) = fillmissing(chloro256_lang_rm(i,:),'previous');
    y(i,:) = nufft(chloro256_lang_rm(i,:),tgrid256(i,:));
end

meanChlorofft = mean(y);

n = length(time200);
f = (0:n-1)/n;

ax14 = figure;
plot(f,abs(y),'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hold on
plot(f,abs(meanChlorofft),'Color','red','DisplayName','Mean Chloropigment Across Depth');
hold off
legend();
xlabel('Normalised Frequency');
ylabel('Power Density (N^4) [s^{-4}]');
title('FFT: chl-a(p,t) at HOTS (1989-2021)');
exportgraphics(ax14,'figures/chla-1988-2021_fft.png');
clear ax14;

%% Simple Gaussian Model of DCM: Platt et al (1988)
% Equation 2.1 from B. Barone (2009), based on Platt et al. (1988)

C0 = 0.01; h = 45; zm = -105; sigma_ = 25; z = 0:-5:-200;

Cz = C0 + (h/(sigma_*sqrt(2*pi)))*exp(-(z-zm).^2/2*sigma_.^2);

%% DCM Model: Uitz et al (2006)
% This adds a slope to the background concentration to the model above.

s = -0.0003;

CzU = C0 + s*z + (h/(sigma_*sqrt(2*pi)))*exp(-(z-zm).^2/2*sigma_.^2);

%% DCM Model Plot

% figure;
% plot(Cz,z,'DisplayName','Platt et al. (1988)');
% hold on
% plot(CzU,z,'DisplayName','Uitz et al (2006)');
% plot(chloro(1:101),-pressure(1:101),'DisplayName','Day One CTD');
% plot(chloro(102:202),-pressure(102:202),'DisplayName','Day Two CTD');
% hold off
% legend();