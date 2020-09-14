function [rotCrystalMatrix] = CrystalMatrixRotZd(crystalMatrix, theta)
%CrystalMatrixRotZd() rotates the crystal matrix about z axis by theta in
%degree.
% Input:
%   crystalMatrix -- crystal matrix, format: [type; proportion; X; Y; Z];
%   theta -- rotation angle in radian;
% Output:
%   rotCrystalMatrix -- rotated crystal matrix;

rotMat = [1, 0, 0, 0, 0;
    0, 1, 0, 0, 0;
    0, 0, cosd(theta), -sind(theta), 0;
    0, 0, sind(theta), cosd(theta), 0;
    0, 0, 0, 0, 1];

rotCrystalMatrix = rotMat * crystalMatrix;

end

