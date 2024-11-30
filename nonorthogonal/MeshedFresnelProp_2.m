function [outWave] = MeshedFresnelProp_2(inWave, fxMesh, fyMesh, wavLen, ...
    propDist, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

propKernel = MeshedFresnelPropKernel_2(fxMesh, fyMesh, wavLen, propDist, ...
    varargin{:});
outWave = ifftshift(ifft2(fftshift(propKernel) .* fft2(fftshift(inWave))));

end