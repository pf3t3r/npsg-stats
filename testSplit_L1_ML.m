clear; clc; close all; addpath("baroneRoutines\");

p_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,4);
chl_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,5);
id_hplc = importdata('data/hplcChla_88-21_110.txt').data(:,1);

mld = importdata("datafiles\MLD.mat").MLD;

%% set up max MLD

testMe = importdata('datafiles\ctd_iso_ALL.mat').ctd;
testMe2 = nan(329,1);
for i = 1:329
    if ~isnan([testMe(i).mld003])
    testMe2(i) = max([testMe(i).mld003]);
    end
end

maxMldPerCruise = testMe2;

%% cruise 1

% so CRN1 MLD = 42 dbar
tmp = num2str(id_hplc);
bottleCRN = str2num(tmp(:,1:3));

testPnew = nan(1,2090);
testCRNnew = nan(1,2090);
testChlNew = nan(1,2090);
testId = nan(1,2090);

testPnew2 = nan(1,2090);
testCRNnew2 = nan(1,2090);
testChlNew2 = nan(1,2090);
testId2 = nan(1,2090);
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
    % with MAX MLD per cruise
    tmpMld2 = maxMldPerCruise(bottleCRN(i));
    if p_hplc(i) < tmpMld2
        testPnew2(i) = p_hplc(i);
        testCRNnew2(i) = bottleCRN(i);
        testChlNew2(i) = chl_hplc(i);
        testId2(i) = id_hplc(i);
%         continue...
    end
end

newP = testPnew(~isnan(testPnew));
newCrn = testCRNnew(~isnan(testCRNnew));
newChl = testChlNew(~isnan(testChlNew));
newId = testId(~isnan(testId));

newChlVec = [newCrn; newId; newP; newChl]';

newP2 = testPnew2(~isnan(testPnew2));
newCrn2 = testCRNnew2(~isnan(testCRNnew2));
newChl2 = testChlNew2(~isnan(testChlNew2));
newId2 = testId2(~isnan(testId2));

newChlVec2 = [newCrn2; newId2; newP2; newChl2]';

%%
figure;
scatter(newChl,newP);

figure;
scatter(newChl2,newP2);
%% clean, bin

% [pb5_hplc,pb10_hplc,chlOut_hplc,n5_hplc,n10_hplc] = cleanAndBin(p_hplc,chl_hplc,id_hplc);       % Chlorophyll a (HPLC method)
[pb5_hplc,pb10_hplc,chlOut_hplc,n5_hplc,n10_hplc] = cleanAndBin(newP,newChl,newId');       % Chlorophyll a (HPLC method)

[pb5_hplc2,pb10_hplc2,chlOut_hplc2,n5_hplc2,n10_hplc2] = cleanAndBin(newP2,newChl2,newId2');       % Chlorophyll a (HPLC method)

%% 10 dbar bin KS test

[ksHp10, obsHp10, dHp10, skHp10, kuHp10, sdHp10, c95hp10, muHp10] = ksOfBinnedCon(chlOut_hplc,pb10_hplc,10);  % 10 dbar / HPLC Chlorophyll a

[ksHp102, obsHp102, dHp102, skHp102, kuHp102, sdHp102, c95hp102, muHp102] = ksOfBinnedCon(chlOut_hplc2,pb10_hplc2,10);  % 10 dbar / HPLC Chlorophyll a

%%

ax1 = figure;
plotKs(dHp10,ksHp10,obsHp10,skHp10,kuHp10,0.5,20.5,true);
sgtitle('HPLC Method: [chl a] (Eulerian Bottle, 10 dbar bin)');
exportgraphics(ax1,'figures/ks_HplcBottleEulerian10db_meanMLD.png'); clear ax1;

%%
ax2 = figure;
plotKs(dHp102,ksHp102,obsHp102,skHp102,kuHp102,0.5,20.5,true);
sgtitle('HPLC Method: [chl a] (Eulerian Bottle, 10 dbar bin)');
exportgraphics(ax2,'figures/ks_HplcBottleEulerian10db_maxMLD.png'); clear ax2;

%% save things

pMldId = newId2;
pMld = newP2;

save testPMld.mat pMld pMldId maxMldPerCruise;