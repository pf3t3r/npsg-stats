clear;
% HPLC Chl-a
% tmp = importdata("data/L1/hplcChla_88-21_150.txt");
% tmp = importdata("data/L0/pc_89-22_200.txt");
tmp = importdata("data/L0/hplcChla_88-22_200.txt");

UT = tmp.data(:,3);

UT(UT<0) = nan;

for i = 1:length(UT)
    HT(i) = UT(i) - 100000;
end

for i = 1:length(UT)
    if HT(i) < 0
        tmp = 240000 + HT(i);
        HT(i) = tmp;
    end
end

time = [UT HT'];

figure;
histogram(HT');

dataInDay = 121/length(UT);