function [grating] = ApertureGrating_X(lx, ly, nx, ny, wavLen, p, p0, theta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

fx = InitFreqAxis(lx, nx);
fy = InitFreqAxis(ly, ny);
[fxMesh, fyMesh] = meshgrid(fx, fy);
rotFxMesh = fxMesh * cosd(theta) - fyMesh * sind(theta);
rotFxMesh = rotFxMesh * wavLen * 1e3 - p0;

grating = (mod(rotFxMesh, p) < p / 2);

end