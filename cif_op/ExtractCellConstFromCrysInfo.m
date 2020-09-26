%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2020  Francis Black Lee and Li Xian

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
function [cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo)
%ExtractCellConstFromCrysInfo() extracts cell constants from the crysInfo 
%struct output by LoadCif().
% Input:
%   crysInfo -- crysInfo struct output by LoadCif();
% Output:
%   cellLengths -- [a, b, c];
%   cellAngles -- [alpha, beta, gamma];

cellLengths = zeros(1, 3);
cellAngles = zeros(1, 3);

% cell length a
cellLengthAIdx = find(strcmp('_cell_length_a', crysInfo.valuedProperty(:, 1)));
cellLengths(1) = str2double(crysInfo.valuedProperty{cellLengthAIdx, 2});

% cell length b
cellLengthBIdx = find(strcmp('_cell_length_b', crysInfo.valuedProperty(:, 1)));
cellLengths(2) = str2double(crysInfo.valuedProperty{cellLengthBIdx, 2});

% cell length c
cellLengthCIdx = find(strcmp('_cell_length_c', crysInfo.valuedProperty(:, 1)));
cellLengths(3) = str2double(crysInfo.valuedProperty{cellLengthCIdx, 2});

% cell angle alpha
cellAngleAlphaIdx = find(strcmp('_cell_angle_alpha', crysInfo.valuedProperty(:, 1)));
cellAngles(1) = str2double(crysInfo.valuedProperty{cellAngleAlphaIdx, 2});

% cell angle beta
cellAngleBetaIdx = find(strcmp('_cell_angle_beta', crysInfo.valuedProperty(:, 1)));
cellAngles(2) = str2double(crysInfo.valuedProperty{cellAngleBetaIdx, 2});

% cell angle gamma
cellAngleGammaIdx = find(strcmp('_cell_angle_gamma', crysInfo.valuedProperty(:, 1)));
cellAngles(3) = str2double(crysInfo.valuedProperty{cellAngleGammaIdx, 2});

end

