function [rotCrystalMatrix] = CrystalMatrixRotX(crystalMatrix, theta)
%CrystalMatrixRotX() rotates the crystal matrix about x axis by theta in 
%radian.
% Input:
%   crystalMatrix -- crystal matrix, format: [type; proportion; X; Y; Z];
%   theta -- rotation angle in radian;
% Output:
%   rotCrystalMatrix -- rotated crystal matrix;

rotMat = [1, 0, 0, 0, 0;
    0, 1, 0, 0, 0;
    0, 0, 1, 0, 0;
    0, 0, 0, cos(theta), -sin(theta);
    0, 0, 0, sin(theta), cos(theta)];

rotCrystalMatrix = rotMat * crystalMatrix;

end

