close all; clc; clear;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 29 15]);
set(0,'defaultAxesFontSize',12);
addpath("baroneRoutines\");

%% Load data

f_e = load('datafiles\f_all.mat').f_er;
f_l = load('datafiles\f_all.mat').f_lag;
t = load('datafiles\f_all.mat').t_copy;

%% Test sunrise/sunset routine
% doesn't work properly... trying routine that BB sent.

lat = 22.75; lon = -158; 
% date = datetime('today');
date = t;

%%
rs = nan(4347,2);

for i = 1:4347
    [tmp,~,~,~,~,~] = suncycle(lat,lon,date(i));

    tmp = tmp - 10;
    for j = 1:2
        if tmp(j) < 0
            tmp(j) = tmp(j) + 24;
        end
    end
    rs(i,:) = tmp;
end

%%
rs2 = hours(rs);
rs2.Format = 'hh:mm';

%% 

% print 1 if that cast is taken at night time

sunriseTime = hours(rs2(:,1)');
sunsetTime = hours(rs2(:,2)');
% test = datenum(t);
dayFrac = rem(datenum(t),1);
castTime = dayFrac*24;

castAtNight = nan(4347,1);

for i = 1:4347
    if castTime(i) < sunriseTime(i) || castTime(i) > sunsetTime(i)
        castAtNight(i) = 1;
    end
end

%%
