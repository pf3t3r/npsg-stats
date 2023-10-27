clc; close all; clear;
addpath("C:\Users\pfarrell\OneDrive - NIOZ\Documenten\hotFTPupcast");
%% Extract DCM from upcast

for i = 1:31
    [pcm(i),sigcm(i),~] = hotFTPupcast(300,i);
end

%% Compare with hotFTPextract dcm

dcmE = load("dcm.mat").dcm;

%% CRN 300
a = 9270; b = 9301;

ax = figure;
plot(dcmE(a:b,2),dcmE(a:b,3),DisplayName='EXTRACT dcm');
hold on
plot(1:1:31,pcm,DisplayName='UPCAST dcm');
hold off
legend(); xlabel('Cast'); ylabel('Pressure [dbar]');
set(gca,"YDir","reverse");
title('HOT-300: Comparison of DCM Calculation Routines');
exportgraphics(ax,'figures/dcmComp/hot300.png');

%% Output no. of casts for each cruise
% There are 31 - lol

noOfCasts = [];
cast = dcmE(:,2);
crn = dcmE(:,1);

tmp = 0;
% tmp2 = nan(3,1);
for i = 1:329
    for j = 1:10199
        if crn(j) == i && ~isnan(cast(j))
            tmp = tmp + 1;
        end
        tmp2(i) = tmp;
    end
end

tmp3 = nan(329,1);
tmp3(1) = 17;
for k = 2:329
    tmp3(k) = tmp2(k) - tmp2(k-1);
end

% tmp3 = no. of casts per cruise

%% find upcast dcms for all cruises

% dcmOutputU = nan(329*31,1);
% for i = 1:329
%     tmpCastNo = tmp3(i);
%     for j = 1:31
%         
%     end
% end