clear; clc; close all; addpath("baroneRoutines\");

p_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,4);
chl_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,5);
id_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,1);

mld = importdata("datafiles\MLD.mat").MLD;

%% cruise 1

% so CRN1 MLD = 42 dbar
tmp = num2str(id_hplc);
bottleCRN = str2num(tmp(:,1:3));

testPnew = nan(1,2090);
testCRNnew = nan(1,2090);
testChlNew = nan(1,2090);
testId = nan(1,2090);
% tmpMld = [];

for i = 1:2062
    %bottleCRN(i)
    tmpMld = mld(bottleCRN(i));
    disp(tmpMld);
    if p_hplc(i) < 2*tmpMld - 2
        testPnew(i) = p_hplc(i);
        testCRNnew(i) = bottleCRN(i);
        testChlNew(i) = chl_hplc(i);
        testId(i) = id_hplc(i);
    end
end

newP = testPnew(~isnan(testPnew));
newCrn = testCRNnew(~isnan(testCRNnew));
newChl = testChlNew(~isnan(testChlNew));
newId = testId(~isnan(testId));

newChlVec = [newCrn; newId; newP; newChl]';
%%
figure;
scatter(newChl,newP);
%% clean, bin

% [pb5_hplc,pb10_hplc,chlOut_hplc,n5_hplc,n10_hplc] = cleanAndBin(p_hplc,chl_hplc,id_hplc);       % Chlorophyll a (HPLC method)
[pb5_hplc,pb10_hplc,chlOut_hplc,n5_hplc,n10_hplc] = cleanAndBin(newP,newChl,newId');       % Chlorophyll a (HPLC method)


%% 10 dbar bin KS test

[ksHp10, obsHp10, dHp10, skHp10, kuHp10, sdHp10, c95hp10, muHp10] = ksOfBinnedCon(chlOut_hplc,pb10_hplc,10);  % 10 dbar / HPLC Chlorophyll a

%%

ax1 = figure;
plotKs(dHp10,ksHp10,obsHp10,skHp10,kuHp10,0.5,20.5,true);
sgtitle('HPLC Method: [chl a] (Eulerian Bottle, 10 dbar bin)');
exportgraphics(ax1,'figures/ks_HplcBottleEulerian10db.png'); clear ax1;