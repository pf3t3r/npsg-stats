function [ks, obs, depth2, Sk, ku, sd, rV, pV] = ksOfBinnedCon(X, p, binning, threshold)
%ksOfBinnedCon find the KS statistic
% INPUTS:
% X = substance concentration,
% p = binned pressure [dbar]
% binning = range of depths to bin values (default=10)
% threshold = no. of values needed for us to consider results (default=100)
% OUTPUTS: 
% ks = KS test for five distributions 
% obs = observations per depth
% depth2 = array of depths above threshold (=100)
% Sk = skewness
% ku = kurtosis
% sd = standard deviation (3 x n). 1: STD of mle (norm); 2: STD of data
% (norm); 3: stdMle/stdData (norm)

if nargin <4
    threshold = 50;
end

if nargin <3
    binning = 10;
end

if binning == 5
    ks = nan(5,40); obs = nan(40,1); n = 40; depth = 5:5:200;
elseif binning == 10
    ks = nan(5,20); obs = nan(20,1); n = 20; depth = 5:10:200;
else
    msg = 'Binning input not valid. Must be either 5 or 10 dbar.';
    error(msg);
end

std = nan(1,n);
for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = X(p==i);
    % apply KS test to X_i
    % change limit below to >3 to fix error with picoeu -> may change other
    % results
    if length(X_i) > 3
        %test = std(X_i);
        %disp(i);
        [~,ks(:,i),~,~,sd(i,:),~] = statsplot2(X_i,'noplot');
        [rV(:,i),pV(:,i)] = bbvuong(X_i);
        %[~,ks(:,i),~,~,~,~] = statsplot2(X_i,'noplot');
        %disp(size(tmpC95));
        Sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
%         tmpDat = [std(X_i) std(log(X_i))];
%         tmpDatMu = [mean(X_i) mean(log(X_i))];
%         tmpComp = [tmpMle(1)/tmpDat(1) tmpMle(2)/tmpDat(2)];
%         tmpCompMu = [muMle(1)/tmpDatMu(1) muMle(2)/tmpDatMu(2)];
%         mu(i,:) = [muMle(1) muMle(2) tmpDatMu(1) tmpDatMu(2) tmpCompMu(1) tmpCompMu(2)];
%         sd(i,:) = [tmpMle(1) tmpMle(2) tmpDat(1) tmpDat(2) tmpComp(1) tmpComp(2)];
%         c95(i,1:2) = [tmpC95(1,1) tmpC95(2,1)];
%         c95(i,3:4) = [tmpC95(1,2) tmpC95(2,2)];
%         c95(i,5:6) = [tmpC95(1,3) tmpC95(2,3)];
%         c95(i,7:8) = [tmpC95(1,4) tmpC95(2,4)];
    end
    obs(i) = length(X_i);
    clear X_i;
end

for i = 1:n
    if obs(i) < threshold
        ks(:,i) = nan;   
        Sk(i) = nan;
        ku(i) = nan;
        sd(i,:) = nan;
        rV(:,i) = nan;
        pV(:,i) = nan;
%         sd(i,:) = nan;
%         c95(i,:) = nan;
%         mu(i,:) = nan;
    end
end

tmp = [];
for i = 1:n
    if ~isnan(sum(ks(:,i)))
        tmp = [tmp i];
    end
end
depth2 = depth(tmp);
Sk = Sk(tmp);
ku = ku(tmp);
sd = sd(tmp,:);
rV = rV(:,tmp);
pV = pV(:,tmp);
% sd = sd(tmp,:);
% c95 = c95(tmp,:);
% mu = mu(tmp,:);

ks = ks(:,~all(isnan(ks)));

end