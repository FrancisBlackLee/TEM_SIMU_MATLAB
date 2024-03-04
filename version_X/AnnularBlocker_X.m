function [blocker] = AnnularBlocker_X(lx, ly, nx, ny, x0, y0, r0, r1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = InitAxis(lx, nx);
y = InitAxis(ly, ny);
[xMesh, yMesh] = meshgrid(x, y);
rMesh = sqrt((xMesh - x0).^2 + (yMesh - y0).^2);
blocker = (rMesh > r0 & rMesh < r1);

end