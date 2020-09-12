function [rotCrystalMatrix] = CrystalMatrixRotYd(crystalMatrix, theta)
%CrystalMatrixRotYd() rotates the crystal matrix about y axis by theta in
%degree.
% Input:
%   crystalMatrix -- crystal matrix, format: [type; proportion; X; Y; Z];
%   theta -- rotation angle in degree;
% Output:
%   rotCrystalMatrix -- rotated crystal matrix;

rotMat = [1, 0, 0, 0, 0;
    0, 1, 0, 0, 0;
    0, 0, cosd(theta), 0, sind(theta);
    0, 0, 0, 1, 0;
    0, 0, -sind(theta), 0, cosd(theta)];

rotCrystalMatrix = rotMat * crystalMatrix;

end

