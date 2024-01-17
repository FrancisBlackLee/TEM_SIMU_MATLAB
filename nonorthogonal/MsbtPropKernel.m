function [propKernel] = MsbtPropKernel(fxMesh, fyMesh, ku, wavLen, propDist)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

propKernel = exp(-1i * pi * wavLen * propDist * (fxMesh.^2 + fyMesh.^2) +...
    1i * 2 * pi * propDist * (fxMesh * ku(1) + fyMesh * ku(2)));

end