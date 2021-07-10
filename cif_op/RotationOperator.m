function [rotMat] = RotationOperator(direction)
%RotationOperator() computes the martix form of the rotation operator given
%the view or projection direction.
% Input:
%   direction -- view or projection direction, vector with 3 elements;
% Output:
%   rotMat -- rotation operator matrix for cartesian coordinates:
%       usage: [xp; yp; zp] = rotMat * [x; y; z];

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

newVecZ = reshape(direction, 1, []);
newVecX = [0, newVecZ(3), -newVecZ(2)];
if ~any(newVecX)
    newVecX = [0, 1, 0];
    newVecZ = newVecZ / norm(newVecZ);
else
    % normalization:
    newVecZ = newVecZ / norm(newVecZ);
    newVecX = newVecX / norm(newVecX);
end

newVecY = cross(newVecZ, newVecX);

rotMat = [newVecX; newVecY; newVecZ];

end

