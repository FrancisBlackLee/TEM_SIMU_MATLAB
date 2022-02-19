function [convMat] = ConversionMatrix(cellLengths, cellAngles)
%ConversionMatrix.m computes the conversion matrix of a unit cell.
% Input:
%   cellLengths -- element 1 for cell length a, 2 for cell length b and 3
%       for cell length c;
%   cellAngles -- element 1 for cell angle alpha (between bases a and c);
%       2 for cell angle beta (between bases b and c) and
%       3 for cell angle gamma (between bases a and b);
% Output:
%   convMat -- conversion matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2022  Francis Black Lee (Li Xian)

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

if all(cellLengths) && all(cellAngles)
    a = cellLengths(1);
    b = cellLengths(2);
    c = cellLengths(3);

    alpha = cellAngles(1);
    beta = cellAngles(2);
    gamma = cellAngles(3);

    c1 = c * cosd(beta);
    c2 = c * (cosd(alpha) - cosd(beta) * cosd(gamma)) / sind(gamma);
    c3 = sqrt(c^2 - c1^2 - c2^2);

    convMat = [a, b * cosd(gamma), c1;
        0, b * sind(gamma), c2;
        0, 0, c3];
else
    error('Invalid input');
end

end


