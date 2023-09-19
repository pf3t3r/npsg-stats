function [ks, lillie_p] = quickMleAndKs(input,dist)
%quickMleAndKs Quickly estimate the maximum likelihood (MLE) that a
%distribution fits a given dataset, and calculate a p-value from the
%Kolmogorov-Smirnov test which tells us how well it fits (in this case a
%higher p-value indicates a better fit).

if nargin < 2
    dist = 'logn';
end

for i = 1:5
    MLE(i,:) = mle(input(:,i),'distribution',dist);
    xCdf(i,:) = linspace(min(input(:,i)) - 2*std(input(:,i)), max(input(:,i)) + 2*std(input(:,i)),length(input(:,i)));
    yCdf(i,:) = cdf(dist, xCdf(i,:), MLE(i,1), MLE(i,2));
    [~,ks(i)] = kstest(input(:,i),[xCdf(i,:)' yCdf(i,:)']);
    [~,lillie_p(i)] = lillietest(input(:,i));
end

end