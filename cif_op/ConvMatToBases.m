function [bases] = ConvMatToBases(convMat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

bases.a = convMat(:, 1)';
bases.b = convMat(:, 2)';
bases.c = convMat(:, 3)';

end

