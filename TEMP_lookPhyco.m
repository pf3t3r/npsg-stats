clc;close all;clear;

x = importdata('data\L1\pe10_00-08_150.txt').data;

phyco = x(:,5); p = x(:,4); t = x(:,2);

figure
scatter(phyco,p,[],"red",'.'); set(gca,"YDir","reverse");
xlabel('Phycoerythrin (10u) [ng/l]'); ylabel('Pressure [dbar]');

ticker = [];
for i = 1:880
    if t(i+1) > t(i)
        ticker = [ticker i+1];
    end
end
test = max(diff(ticker));
phycoA = nan(test,length(ticker));
tA = nan(size(phycoA));
pA = nan(size(phycoA));


% j = 1
tA(1:11,1) = t(1:11); pA(1:11,1) = p(1:11); phycoA(1:11,1) = phyco(1:11);
for j = 2:75
    disp(j);
    tA(1:(ticker(j)-ticker(j-1)),j) = t(ticker(j-1):ticker(j)-1);
    pA(1:(ticker(j)-ticker(j-1)),j) = p(ticker(j-1):ticker(j)-1);
    phycoA(1:(ticker(j)-ticker(j-1)),j) = phyco(ticker(j-1):ticker(j)-1);
end

% tA2 = datetime(tA, 'ConvertFrom', 'mmddyy', 'Format','MM/dd/yyyy');
%%

figure
contourf(tA(:,1:6),pA(:,1:6),phycoA(:,1:6),linspace(0,101.13,100),"LineColor","auto"); set(gca,"YDir","reverse");
colormap(flipud(cbrewer2('Spectral',100)));

%%
d = string(tA);
for i = 1:75
    for j = 1:22
        if strlength(d(j,i)) < 6
            disp(i);
            d(:,i) = "0" + d(:,i);
        end
    end
end