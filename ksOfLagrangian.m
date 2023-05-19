function [trange,ks,obsPerBin,Sk,Ku,c95] = ksOfLagrangian(id,p,dcmArray,X,Ltid,threshold)
%funShit3 quickly find the DCM-centred (Lagrangian) transformation for a
%given variable.
% INPUTS:
% id: bottle ID
% p = pressure
% dcmArray = shows where the DCM is
% X = bottle concentration
% Ltid = not sure how this works but hey it works...
% OUTPUTS:
% bottleArray = not needed anymore, to remove
% trange = depths where sufficient measurements are present
% chl = also not used??
% pB = binned pressure, also not used?
% ks = Kolmogorov-Smirnov Statistic (p value). High p-value indicates that
% the given distribution fits the data better.
% obsPerBin = no. of observations in a particular depth bin. 


if nargin < 6
    threshold = 100;
end

CRN = str2num(id(:,1:3)); 
cast = str2num(id(:,6:8));
cast(cast==100) = nan;

bottleArray = [CRN cast p];

t = rmmissing(unique(bottleArray(:,1:2),"rows"));

dcmArrayRowNo = [];

for i = 1:length(dcmArray(:,1))
    for x = 1:length(t)
        if dcmArray(i,1:2) == t(x,1:2) 
            dcmArrayRowNo = [dcmArrayRowNo i];
        end
    end
end

% save when bottleArray changes
tid = [];
for i = 2:length(p)
    if bottleArray(i,1) > bottleArray(i-1,1) || bottleArray(i,2) > bottleArray(i-1,2)
        tid = [tid i];
    end
end

tPcm = nan(length(p),1);
tPcm(1:tid(1)-1) = dcmArray(dcmArrayRowNo(1),3);
tPcm(tid(end):end) = dcmArray(dcmArrayRowNo(end),3);

for i = 2:Ltid-2
    disp(i);
    tPcm(tid(i):tid(i+1)-1) = dcmArray(dcmArrayRowNo(i),3);
end

bottleArray = [bottleArray tPcm];

tPLagrangian = nan(length(p),1);
tPLagrangian = bottleArray(:,3) - bottleArray(:,4);
bottleArray = [bottleArray tPLagrangian];

pB10 = round(bottleArray(:,5),-1);
bottleArray = [bottleArray pB10];

bottleArray = [bottleArray X];
tmin = min(bottleArray(:,6));
tmax = max(bottleArray(:,6));
trange = tmin:10:tmax;

chl = bottleArray(:,7);
pB = bottleArray(:,6);
ks = nan(5,length(trange));
obsPerBin = nan(1,length(trange));

Sk = nan(1,length(trange));

for i = 1:length(trange)
    tmp = chl(pB==trange(i));
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obsPerBin(i) = length(tmp);
    if length(tmp) > 3
        disp(i);
        [~,ks(:,i),~] = statsplot2(tmp,'noplot');
        Sk(i) = skewness(tmp);
        Ku(i) = kurtosis(tmp);
    end
end

for i = 1:length(trange)
    if obsPerBin(i) < threshold
        ks(:,i) = nan;
        Sk(i) = nan;
        Ku(i) = nan;
    end
end

end