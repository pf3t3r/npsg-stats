clear; clc; close all;
addpath("baroneRoutines\");

%% Load Test Data
chloro = load("datafiles\chloro.mat"','chloro256_n').chloro256_n;
time = load("datafiles\chloro.mat","time256").time256;

%% Remove NaN or Zeros

tickyboi = [];
antiTick = [];

for i = 1:329
    for j = 1:129
        if (~any(chloro(j,i)))
            tickyboi = [tickyboi i];
            break;
        else
            antiTick = [antiTick i];
            break;
        end
    end
end

%% Create Arrays of chl-a and time which have no zero or NaN data
chloroNew = chloro(:,[antiTick]);
timeNew = time(:,[antiTick]);

%% Look

R = []; p2 = [];

for i = 1:20
%     disp(i);
    [R(:,i),p2(:,i)] = bbvuong(chloroNew(i,:));
end


%% show

figure
plot(R(1,:),0:-1:-20);

% [R,p2] = bbvuong(chloroNew(55,:));