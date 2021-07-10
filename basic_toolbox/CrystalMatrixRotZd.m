function [rotCrystalMatrix] = CrystalMatrixRotZd(crystalMatrix, theta)
%CrystalMatrixRotZd() rotates the crystal matrix about z axis by theta in
%degree.
% Input:
%   crystalMatrix -- crystal matrix, format: [type; proportion; X; Y; Z];
%   theta -- rotation angle in radian;
% Output:
%   rotCrystalMatrix -- rotated crystal matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rotMat = [1, 0, 0, 0, 0;
    0, 1, 0, 0, 0;
    0, 0, cosd(theta), -sind(theta), 0;
    0, 0, sind(theta), cosd(theta), 0;
    0, 0, 0, 0, 1];

rotCrystalMatrix = rotMat * crystalMatrix;

end

