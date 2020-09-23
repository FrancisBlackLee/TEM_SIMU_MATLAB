function [rotMat] = RotationOperator(direction)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

newVecZ = reshape(direction, 1, []);
newVecX = [0, newVecZ(3), -newVecZ(2)];
% normalization:
newVecZ = newVecZ / norm(newVecZ);
newVecX = newVecX / norm(newVecX);

newVecY = cross(newVecZ, newVecX);

rotMat = [newVecX; newVecY; newVecZ];

end

