function [idOut,pOut,XOut] = extractSMLC(id,p,X,pMaxMld)
% This function ...
% INPUTS
% id
% p
% X
% pMaxMld
% OUTPUTS
% idOut
% pOut
% XOut

X(X==-9) = nan;

id = id(~isnan(X),:);
p = p(~isnan(X));
X1 = X(~isnan(X));

crn = str2num(id(:,1:3));

for i = 1:length(crn)
    if crn(i) == 330
        stop = i;
        break
    else
        stop = length(p) + 1;
    end
end

% Extract meas below pMaxMld
L = stop-1; % No. of casts with chl measurements in cruise 1 - 329

tmpP_subML = nan(L,1);
tmpCRN_subML = nan(L,1);
tmpX_subML = nan(L,1);
botID = [];

% To find the measurements taken beneath pMaxMld, we save them to a new
% file...
for i = 1:L
    %crn(i)
    % with MAX MLD per cruise
    tmpMld = pMaxMld(crn(i));
    if p(i) > tmpMld
        tmpP_subML(i) = p(i);
        tmpCRN_subML(i) = crn(i);
        tmpX_subML(i) = X1(i);
        tmpStr = id(i,:);
        botID = [botID;tmpStr];
    end
end

% ...and remove nan values (which represent measurements above pMaxMld)
pOut = tmpP_subML(~isnan(tmpP_subML));
XOut = tmpX_subML(~isnan(tmpX_subML));
idOut = botID;

end