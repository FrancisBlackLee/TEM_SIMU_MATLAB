function [basis] = CrystalIndicesToBasis(convMat, indices)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

basis = indices(1) * convMat(:, 1) +...
        indices(2) * convMat(:, 2) +...
        indices(3) * convMat(:, 3);

end


