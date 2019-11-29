%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

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
function [SliceOut] = SquareLattExpan(SliceIn, LattConst, CellNum)
%SquareLattExpan.m expands the sliced unit cell lattice by the input
%CellNum [Nx, Ny] and output the coordinates. Notice that the input SliceIn
%is the proportional coordinates of atoms in a unit cell or basic square
%lattice.
%   SliceIn -- proportional atomic coordinates (scaled to 0~1) [x1,..., xN;
%               y1,..., yN];
%   LattConst -- planar lattice constants [a, b];
%   SliceDist -- distances between each two slices [D1,..., DN];
%   CellNum -- expansion numbers by which the SliceIn is expanded [Nx, Ny];

% Transpose the matrix to improve speed
SliceOut = SliceIn';
CoordShift = CellNum / 2 + 1; % center the lattice
DistError = 1e-2;
% Expand along x
SliceBase = SliceOut;
for i = 1: CellNum(1) + 1
    SliceOut = [SliceOut; [SliceBase(:, 1) + i, SliceBase(:, 2)]];
end
SliceBase = SliceOut;
for i = 1: CellNum(2) + 1
    SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2) + i]];
end
SliceOut(:, 1) = SliceOut(:, 1) - CoordShift(1);
SliceOut(:, 2) = SliceOut(:, 2) - CoordShift(2);
[row, column] = find((SliceOut(:, 1) <= CellNum(1) / 2 + DistError) & (SliceOut(:, 1) >= -CellNum(1) / 2 - DistError) ...
                   & (SliceOut(:, 2) <= CellNum(2) / 2 + DistError) & (SliceOut(:, 2) >= -CellNum(2) / 2 - DistError));
SliceOut = SliceOut(row, :);
SliceOut = uniquetol(SliceOut, DistError, 'ByRows', true);
SliceOut(:, 1) = LattConst(1) * SliceOut(:, 1);
SliceOut(:, 2) = LattConst(2) * SliceOut(:, 2);
SliceOut = SliceOut';

end

