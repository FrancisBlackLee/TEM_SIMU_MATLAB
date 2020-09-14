function [rotCrystalMatrix] = CrystalMatrixRotZ(crystalMatrix, theta)
%CrystalMatrixRotZ() rotates the crystal matrix about z axis by theta in
%radian.
% Input:
%   crystalMatrix -- crystal matrix, format: [type; proportion; X; Y; Z];
%   theta -- rotation angle in radian;
% Output:
%   rotCrystalMatrix -- rotated crystal matrix;

rotMat = [1, 0, 0, 0, 0;
    0, 1, 0, 0, 0;
    0, 0, cos(theta), -sin(theta), 0;
    0, 0, sin(theta), cos(theta), 0;
    0, 0, 0, 0, 1];

rotCrystalMatrix = rotMat * crystalMatrix;

end

