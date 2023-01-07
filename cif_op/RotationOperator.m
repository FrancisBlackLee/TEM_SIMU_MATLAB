function [rotMat] = RotationOperator(directionZ, directionX)
%RotationOperator() computes the martix form of the rotation operator given
%the view or projection direction.
% Input:
%   direction -- view or projection direction, vector with 3 elements;
% Output:
%   rotMat -- rotation operator matrix for cartesian coordinates:
%       usage: [xp; yp; zp] = rotMat * [x; y; z];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

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

newVecZ = reshape(directionZ, 1, []);
newVecZ = newVecZ / norm(newVecZ);

if nargin == 1
    newVecX = [0, newVecZ(3), -newVecZ(2)];
    if ~any(newVecX)
        newVecX = [0, 1, 0];
        newVecY = [0, 0, 1];
    else
        % normalization:
        newVecX = newVecX / norm(newVecX);
        newVecY = cross(newVecZ, newVecX);
    end
elseif nargin == 2
    newVecX = reshape(directionX, 1, []);
    tolerance = 1.0e-8;
    if abs(dot(newVecX, newVecZ)) < tolerance
        newVecX = newVecX / norm(newVecX);
        newVecY = cross(newVecZ, newVecX);
    else
        error('Input direction x is not vertical to direction z');
    end
else
    error('Invalid input');
end

rotMat = [newVecX; newVecY; newVecZ];

end



