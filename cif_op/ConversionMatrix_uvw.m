function [convMat] = ConversionMatrix_uvw(cellLengths, cellAngles, uvw, hkl)
%ConversionMatrix_uvw() computes the conversion matrix, given the cell
%constants and zone axis / orientation indices.
% Input:
%   cellLengths -- element 1 for cell length a, 2 for cell length b and 3
%       for cell length c;
%   cellAngles -- element 1 for cell angle alpha (between bases b and c);
%       2 for cell angle beta (between bases a and c) and
%       3 for cell angle gamma (between bases a and b);
%   uvw -- zone axis / orientation indices;
% Output:
%   convMat -- conversion matrix;

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

if all(cellLengths) && all(cellAngles) && any(uvw)
    initConvMat = ConversionMatrix(cellLengths, cellAngles);

    viewDirection = CrystalIndicesToBasis(initConvMat, uvw);

    if nargin == 3
        rotMat = RotationOperator(viewDirection);
    elseif nargin == 4
        horizonalDirection = CrystalIndicesToBasis(initConvMat, hkl);
        rotMat = RotationOperator(viewDirection, horizonalDirection);
    else
        error('Incorrect number of input arguments.');
    end

    convMat = rotMat * initConvMat;
else
    convMat = zeros(3);
end

end



