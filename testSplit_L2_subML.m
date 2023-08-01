clear; clc; close all; addpath("baroneRoutines\");

p_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,4);
chl_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,5);
id_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,1);

% mld = importdata("datafiles\MLD.mat").MLD;

%%

% Only look at values below the maximum MLD of a cruise
% We don't have a convenient measurement of MLD in 'isopycnal' coordinates.

pMldId = load('testPMld.mat').pMldId;
pMld = load('testPMld.mat').pMld;
pMaxMld = load('testPMld.mat').maxMldPerCruise;


tmp = num2str(id_hplc);
bottleCRN = str2num(tmp(:,1:3));

%%

% these new values are pressures etc that are below MLD

for i = 1:2062
    %bottleCRN(i)
    % with MAX MLD per cruise
    tmpMld = pMaxMld(bottleCRN(i));
    if p_hplc(i) > tmpMld
        testPnew2(i) = p_hplc(i);
        testCRNnew2(i) = bottleCRN(i);
        testChlNew2(i) = chl_hplc(i);
        testId2(i) = id_hplc(i);
%         continue...
    end
end