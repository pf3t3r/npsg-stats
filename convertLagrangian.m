function [lagrangianData] = convertLagrangian(eulerianData,offset,n)
%convertLagrangian converts data in depth/pressure coordinates into
%Lagrangian coordinates, i.e. coordinates that are centered on the DCM.
% INPUT: eulerianData, offset, n = length of dataset
% OUTPUT: lagrangianData

for k = 1:length(n)
    lagrangianData(:,k) = circshift(eulerianData,offset);
    if offset > -1 && offset < 40
        lagrangianData(1:offset,k) = NaN;
    elseif offset == -1
        lagrangianData(end,k) = NaN;
    elseif offset < -1 && offset > -40
        lagrangianData((end+offset):end,k) = NaN;
    elseif abs(offset) > 40
        lagrangianData(:,k) = NaN;
    end
end