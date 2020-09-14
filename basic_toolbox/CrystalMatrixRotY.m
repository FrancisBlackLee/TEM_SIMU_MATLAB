function [rotCrystalMatrix] = CrystalMatrixRotY(crystalMatrix, theta)
%CrystalMatrixRotY() rotates the crystal matrix about y axis by theta in
%radian.
% Input:
%   crystalMatrix -- crystal matrix, format: [type; proportion; X; Y; Z];
%   theta -- rotation angle in radian;
% Output:
%   rotCrystalMatrix -- rotated crystal matrix;

rotMat = [1, 0, 0, 0, 0;
    0, 1, 0, 0, 0;
    0, 0, cos(theta), 0, sin(theta);
    0, 0, 0, 1, 0;
    0, 0, -sin(theta), 0, cos(theta)];

rotCrystalMatrix = rotMat * crystalMatrix;

end

