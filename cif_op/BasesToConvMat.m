function [convMat] = BasesToConvMat(bases)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

vecA = reshape(bases.a, [], 1);
vecB = reshape(bases.b, [], 1);
vecC = reshape(bases.c, [], 1);

convMat = [vecA, vecB, vecC];

end

